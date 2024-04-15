/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.


SCHEME CELL LEVEL 0 TO BE SPLIT, NEXT TO ONE CELL OF LEVEL 1

              X------------X|
             /    b       / |
            /            /  |
            X-----0-----X   |
            |     |  a  |   |
           0|_____|0____|0  X    
            |     |     |  /
            |     |     | /
            X-----0-----X/

X: anchor points, of level 0
0: level 1 points, no anchor points (are anchor points of a cell of level 1)

b: face with 4 anchor points (and 5 points)
a: face with 1 anchor point (the other 3 are of higher level than cell).

\*---------------------------------------------------------------------------*/

#include "hexRefPoliMi.H"

#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "polyAddCell.H"
#include "polyModifyFace.H"
#include "syncTools.H"
#include "faceSet.H"
#include "cellSet.H"
#include "pointSet.H"
#include "OFstream.H"
#include "Time.H"
#include "FaceCellWave.H"
#include "mapDistributePolyMesh.H"
#include "refinementData.H"
#include "refinementDistanceData.H"
#include "degenerateMatcher.H"

#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hexRefPoliMi, 0);

    //- Reduction class. If x and y are not equal assign value.
    template<label value>
    class ifEqEqOp
    {
        public:
        void operator()(label& x, const label y) const
        {
            x = (x == y) ? x : value;
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::hexRefPoliMi::is2D()
{
    bool is2D_= false;

    label nEmpty = 0;

    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (isA<emptyPolyPatch>(mesh_.boundaryMesh()[patchi]) ||
        isA<wedgePolyPatch>(mesh_.boundaryMesh()[patchi])
        )
        {
            nEmpty += mesh_.boundaryMesh()[patchi].size();
        }
    }

    reduce(nEmpty, sumOp<label>());

    if (nEmpty)
    {
        is2D_= true;
    }

    return is2D_;
}


/*
-------- MY EXPLINATION of getAnchorCell --------------

We are trying to get the anchor cell of a point.
If a cell is not refined, we obtained celli itself.
In a refined cell we can have 2 situations:

a) pointi is an old vertix (of the big cell):
	we are going to return the celli of the cell that was added in the 
	corner of this point, belonging to celli (among the 8 into which 
	got split, the one in pointi)
b)pointi is NOT an old vertix (is middle of edge):
	as it does not coincide with any original anchor points, we are going
	to explore the 4 points of facei. One of these has to be a vertex of the
	original big cell. Finding this vertix, we can procede as in a)
-----------------------------------------------------*/

// NEW 2D version
// Check whether pointi is an anchor on celli.
// If it is not check whether any other point on the face is an anchor cell.
Foam::label Foam::hexRefPoliMi::getAnchorCell
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label celli,
    const label facei,
    const label pointi
) const
{
    if(is2D_) // PoliMi
    {
        if (cellAnchorPoints[celli].size())
        {
            // list of integers to store the logic relation between the vertices
            // of each cell
            labelList logicRelation(8, -1);
           
            label masterPatch = -1;
            label slavePatch = -1;
            bool firstDone(false);

            // find master (first appearing) and slave patches
            for(int patchi = 0; patchi < mesh_.boundaryMesh().size(); patchi ++)
            {
                if (!firstDone)
                {
                    if (mesh_.boundaryMesh()[patchi].type()=="empty" ||
                        mesh_.boundaryMesh()[patchi].type()=="wedge" )
                    {
                        masterPatch = patchi;
                        firstDone = true;
                    }
                }
                else
                {
                    if (mesh_.boundaryMesh()[patchi].type()=="empty" ||
                        mesh_.boundaryMesh()[patchi].type()=="wedge" )
                    {
                        slavePatch = patchi;
                        break;
                    }
                }
            }

            // find master and slave faces of celli
            label masterFace = -1;
            label slaveFace = -1;

            labelList cFaces = mesh_.cells()[celli];

            forAll(cFaces, i)
            {
                label fi = cFaces[i];
                if(!mesh_.isInternalFace(fi))
                {
                    label patchfi = mesh_.boundaryMesh().whichPatch(fi);

                    if(patchfi == masterPatch)
                    {
                        masterFace = fi;
                    }
                    else if(patchfi == slavePatch)
                    {
                        slaveFace = fi;
                    }
                }
            }

            const face& mFace = mesh_.faces()[masterFace];
            const face& sFace = mesh_.faces()[slaveFace];

            if(debug)
            {
                Pout << "masterPatch: " << masterPatch << endl;
                Pout << "slavePatch: " << slavePatch << endl;
                Pout << "masterFace: " << masterFace << endl;
                Pout << "slaveFace: " << slaveFace << endl;
            }

            forAll(cellAnchorPoints[celli], anchorPosition)
            {
                // we want to find the projection of the anchor point on the master
                // face. goodFace is true if it's on the master face
                bool goodFace = false;

                label anchorPoint = cellAnchorPoints[celli][anchorPosition];

                if(debug)
                {
                    Pout << "cellAnchorPoints " << cellAnchorPoints[celli] << endl;
                    Pout << "anchorPoint " << anchorPoint << endl;
                    Pout << "anchorPosition " << anchorPosition << endl;
                }

                label isMaster = findIndex(mFace, anchorPoint);
                label isSlave = findIndex(sFace, anchorPoint);

                if(isMaster != -1 && isSlave == -1)
                {
                    // anchor is in master face
                    goodFace = true;
                }
                else if(isMaster == -1 && isSlave != -1)
                {
                    goodFace = false;
                    // anchor is in slave face
                }
                else
                {
                    FatalErrorInFunction
                    << "pointi is not on the master or the slave patches"
                    << endl
                    << abort(FatalError);
                }
                
                if(goodFace)
                {
                    // if it's on the master face, we can store it
                    logicRelation[anchorPosition] = anchorPosition;
                }
                else
                {
                    // if it's on the slave face, we find its projection
                    // into the master face
                    const labelList& pEdges = mesh_.pointEdges()[anchorPoint];
                    
                    forAll(pEdges, ei)
                    {
                        label currentEdge = pEdges[ei];
                        const edge& cEdge = mesh_.edges()[currentEdge];

                        label otherVertex = cEdge.otherVertex(anchorPoint);
                        label isInMaster = findIndex(mFace, otherVertex);

                        if(isInMaster != -1)
                        {
                            label projectedPosition = findIndex(cellAnchorPoints[celli], otherVertex);
                            logicRelation[anchorPosition] = projectedPosition;

                            break;
                        }
                    }
                }
            }

            // Re-scale logicRelation to be between 0-3
            labelList mapping(8,-1);
            label counter = 0;

            forAll(logicRelation, i)
            {
                if (mapping[logicRelation[i]] == -1)
                {
                    mapping[logicRelation[i]] = counter++;
                }
                logicRelation[i] = mapping[logicRelation[i]];
            }

            if(debug)
            {
                Pout << "logicRelation: " << logicRelation << endl;
            }

            // return corresponding cell of cellAddedCells
            // To do so, I have to count the number of 
            // non-projected points placed before the one I want to return

            label inputPosition = findIndex(cellAnchorPoints[celli], pointi);

            if (inputPosition != -1)
            {
                label indexValue = logicRelation[inputPosition];
                return cellAddedCells[celli][indexValue];
            }

            // pointi is not an anchor cell.
            // Maybe we are already a refined face so check all the face
            // vertices.
            const face& f = mesh_.faces()[facei];

            forAll(f, fp)
            {
                label inputPosition = findIndex(cellAnchorPoints[celli], f[fp]);

                if (inputPosition != -1)
                {
                    label indexValue = logicRelation[inputPosition]; 
                    return cellAddedCells[celli][indexValue];
                }
            }

            // Problem.
            dumpCell(celli);
            Perr<< "cell:" << celli << " anchorPoints:" << cellAnchorPoints[celli]
                << endl;

            FatalErrorInFunction
                << "Could not find point " << pointi
                << " in the anchorPoints for cell " << celli << endl
                << "Does your original mesh obey the 2:1 constraint and"
                << " did you use consistentRefinement to make your cells to refine"
                << " obey this constraint as well?"
                << abort(FatalError);

            return -1;
        }
        else
        {
        return celli;
        }
    }

    else

    {
        //Copy original code for 3D case from hexRef8
        if (cellAnchorPoints[celli].size())
        {
            label index = findIndex(cellAnchorPoints[celli], pointi);

            if (index != -1)
            {
                return cellAddedCells[celli][index];
            }


            // pointi is not an anchor cell.
            // Maybe we are already a refined face so check all the face
            // vertices.
            const face& f = mesh_.faces()[facei];

            forAll(f, fp)
            {
                label index = findIndex(cellAnchorPoints[celli], f[fp]);

                if (index != -1)
                {
                    return cellAddedCells[celli][index];
                }
            }

            // Problem.
            dumpCell(celli);
            Perr<< "cell:" << celli << " anchorPoints:" << cellAnchorPoints[celli]
                << endl;

            FatalErrorInFunction
                << "Could not find point " << pointi
                << " in the anchorPoints for cell " << celli << endl
                << "Does your original mesh obey the 2:1 constraint and"
                << " did you use consistentRefinement to make your cells to refine"
                << " obey this constraint as well?"
                << abort(FatalError);

            return -1;
        }
        else
        {
            return celli;
        }
    }
}

Foam::label Foam::hexRefPoliMi::storeMidPointInfo
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& edgeMidPoint,
    const label celli,
    const label facei,
    const bool faceOrder,
    const label edgeMidPointi,
    const label anchorPointi,
    const label faceMidPointi,

    Map<edge>& midPointToAnchors,
    Map<edge>& midPointToFaceMids,
    polyTopoChange& meshMod
) const
{
    bool changed = false;
    bool haveTwoAnchors = false;

    Map<edge>::iterator edgeMidFnd = midPointToAnchors.find(edgeMidPointi);

    if (edgeMidFnd == midPointToAnchors.end())
    {
        midPointToAnchors.insert(edgeMidPointi, edge(anchorPointi, -1));
    }
    else
    {
        edge& e = edgeMidFnd();

        if (anchorPointi != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = anchorPointi;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoAnchors = true;
        }
    }

    bool haveTwoFaceMids = false;

    Map<edge>::iterator faceMidFnd = midPointToFaceMids.find(edgeMidPointi);

    if(!is2D_) // PoliMi
    {

        if (faceMidFnd == midPointToFaceMids.end())
        {
            midPointToFaceMids.insert(edgeMidPointi, edge(faceMidPointi, -1));
        }
        else
        {
        
            edge& e = faceMidFnd();

            if (faceMidPointi != e[0])
            {
                if (e[1] == -1)
                {
                    e[1] = faceMidPointi;
                    changed = true;
                }
            }

            if (e[0] != -1 && e[1] != -1)
            {
                haveTwoFaceMids = true;
            }
        }
    }
    else
    {
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // PoliMi
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // In 2D, we will never have haveToFaceMids = true.
        // We will always have (faceMidPointi, -1)
        // We need to loop for the first empty face without doing nothing
        // only to store the index of faceMidPointi in the Map.
        // Then, when entering from the other face, we will have the 2 values
        // of the faceMids in different location of the map.
        // We will have to look for it and set haveTwoFaceMids = true
    
        if (faceMidFnd == midPointToFaceMids.end())
        {
            labelList tableOfContents = midPointToFaceMids.toc();
            label otherFaceMid = -1;

            forAll(tableOfContents, i)
            {
                label toci = tableOfContents[i];
                Map<edge>::iterator iter = midPointToFaceMids.find(toci);

                edge& oldEdge = iter();//this should be (some faceMidPoint, -1)

                if (oldEdge[0] != faceMidPointi)
                {
                    otherFaceMid = oldEdge[0];
                    break;
                }
            }

            midPointToFaceMids.insert(edgeMidPointi, edge(faceMidPointi, otherFaceMid));

            changed = true;
    
        }

        Map<edge>::iterator checkTwoFaceMids = midPointToFaceMids.find(edgeMidPointi);

        edge& e = checkTwoFaceMids();

        if(e[0] != -1 && e[1] != -1)
        {
            haveTwoFaceMids = true;
        }
    }

    // Now, in 2D, I will have the following modPointToFaceMids:
    // - For the edges on the first face, it's going to be (faceMidi, -1), with the index
    // of the center of this first face. Here, we will never have haveTwoFaceMids = true, so we don't
    // create any face
    // - For the edges on the second face, we are going to have couples like:
    // (faceMidi (of the current/second face), faceMidi (of the other/first face))
    // so something like:
    // 0 (8, -1)
    // 1 (8, -1)
    // 2 (8, -1)
    // 3 (8, -1)
    // 4 (9, 8)
    // 5 (9, 8)
    // 6 (9, 8)
    // 7 (9. 8)
    // With 8 and 9 being the face centers that we need to build the internal faces.
    // Now, we just have to find the corresponding edgeMidPoints, one being the projection of the other
    // That is not that easy to do but I think we can use the code I did for the 2D case the other time
    // or think of a better way, closer to the original 3D one.
    // Remember that once the points have been correctly ordered, the creation of the faces uses the same
    // code as in 3D
    // Plus, with this method, we are always introducing the points in the same order, so we might skip
    // the backwards walk. Have always in mind that we call this function twice for each anchor point
    // on the 2 directions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Check if this call of storeMidPointInfo is the one that completed all
    // the necessary information.

    if (changed && haveTwoAnchors && haveTwoFaceMids)
    {
        //const edge& anchors = midPointToAnchors[edgeMidPointi];
        const edge& faceMids = midPointToFaceMids[edgeMidPointi];

        label otherFaceMidPointi = faceMids.otherVertex(faceMidPointi);

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // PoliMi
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // In 2D I need to find the projection of edgeMidPointi on the other face
        // Maybe I can do all of this in a private function that returns otherFaceMidPoint?
        
        label otherEdgeMidPointi = -1;

        if(is2D_)
        {
            const cell& cFaces = mesh_.cells()[celli];

            label otherFace = -1;

            for(int j = cFaces.size()-1; j >=0; j--)
            {
                label facej = cFaces[j];

                if(!mesh_.isInternalFace(facej))
                {
                    label patchID = mesh_.boundaryMesh().whichPatch(facej);

                    if(facej != facei
                        && 
                        (mesh_.boundaryMesh()[patchID].type()=="empty"
                         || 
                         mesh_.boundaryMesh()[patchID].type()=="wedge")
                        )
                    {
                        otherFace = facej;
                        break;
                    }
                }
            }

//            const edge& anchors = midPointToAnchors[edgeMidPointi];
            // both anchors at both sides of edgeMidPointi

            const face& oFace = mesh_.faces()[otherFace];// 'other' face

            const labelList& pEdges = mesh_.pointEdges()[anchorPointi];
            // All edges with the anchorPointi

            label otherPoint = -1;
            forAll(pEdges, ei)
            {
                label currentEdge = pEdges[ei];
                const edge& cEdge = mesh_.edges()[currentEdge];

                otherPoint = cEdge.otherVertex(anchorPointi);
                label otherPointIsInOtherFace = findIndex(oFace, otherPoint);

                if(otherPointIsInOtherFace != -1)
                {
                    break;
                }
            }

            if(otherPoint == -1)
            {
                dumpCell(mesh_.faceOwner()[facei]);
                if (mesh_.isInternalFace(facei))
                {
                    dumpCell(mesh_.faceNeighbour()[facei]);
                }

                FatalErrorInFunction
                << " Could not find projection of edgeMidPointi in 2D on celli = " << celli 
                << abort(FatalError);
            }

            // Now, we have to find the edgeMidPoint from the otherPoint (projection
            // of the anchorPointi) 

            const labelList& ofEdges = mesh_.faceEdges(otherFace);
            if(faceOrder)
            {
                label otherPosition = findIndex(oFace, otherPoint);
                if(pointLevel_[oFace[oFace.rcIndex(otherPosition)]]
                        <= cellLevel_[celli])
                {
                    label otherEdge = ofEdges[oFace.rcIndex(otherPosition)];
                    otherEdgeMidPointi = edgeMidPoint[otherEdge];
                }
                else
                {
                    label otherEdge = findLevel
                            (
                            otherFace, 
                            oFace,
                            oFace.rcIndex(otherPosition),
                            false, 
                            cellLevel_[celli] +1
                            );
                    otherEdgeMidPointi = oFace[otherEdge];
                }
            }
            else
            {
                label otherPosition = findIndex(oFace, otherPoint);

                if(pointLevel_[oFace[oFace.fcIndex(otherPosition)]]
                        <= cellLevel_[celli])
                {
                    label otherEdge = ofEdges[otherPosition];
                    otherEdgeMidPointi = edgeMidPoint[otherEdge];
                }
                else
                {
                    label otherEdge = findLevel
                            (
                            otherFace,
                            oFace,
                            oFace.fcIndex(otherPosition),
                            false,
                            cellLevel_[celli] +1
                            );
                    otherEdgeMidPointi = oFace[otherEdge];
                }
            }
        }

        if(is2D_ && otherEdgeMidPointi == -1)
        {
            FatalErrorInFunction
            << " not found otherEdgeMidPointi:" << facei
            << abort(FatalError);
        }

        // Here, I already have all the points I need to create the internal faces
        // both in 2D and 3D
        // In 3D, otherEdgeMidPointi = -1;
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Create face consistent with anchorI being the owner.
        // Note that the edges between the edge mid point and the face mids
        // might be marked for splitting. Note that these edge splits cannot
        // be between cellMid and face mids.

        const edge& anchors = midPointToAnchors[edgeMidPointi];
 
        // find the projection of anchorPointi

        label otherAnchorPointi = -1;
        if(is2D_) // PoliMi
        {
            const edge& otherAnchors = midPointToAnchors[otherEdgeMidPointi];

            forAll(otherAnchors, i)
            {
                label isProjection = meshTools::findEdge(mesh_, anchorPointi, otherAnchors[i]);
                if (isProjection != -1)
                {
                    otherAnchorPointi = otherAnchors[i];
                }
            }

            if (otherAnchorPointi == -1)
            {
                //error
                FatalErrorInFunction 
                    << "could not find projection of anchorPointi"
                    << abort(FatalError);
            }
        }


        // anchors on the master face, the 1st entering

        DynamicList<label> newFaceVerts(4);

        if (!is2D_ && faceOrder == (mesh_.faceOwner()[facei] == celli))
            // PoliMi : added if not 2D condition
        {
            newFaceVerts.append(faceMidPointi);

            // Check & insert edge split if any
            insertEdgeSplit
            (
                edgeMidPoint,
                faceMidPointi,  // edge between faceMid and
                edgeMidPointi,  // edgeMid
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointi,
                otherFaceMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(otherFaceMidPointi);
            newFaceVerts.append(cellMidPoint[celli]);
        }
        else if(!is2D_ && faceOrder != (mesh_.faceOwner()[facei] == celli))
            // PoliMi added not 2D condition
        {
            newFaceVerts.append(otherFaceMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                otherFaceMidPointi,
                edgeMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointi,
                faceMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(faceMidPointi);
            newFaceVerts.append(cellMidPoint[celli]);
        }
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // PoliMi
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        else if(is2D_ && (faceOrder == (mesh_.faceOwner()[facei] == celli)))
        {
            newFaceVerts.append(faceMidPointi);

            // Check & insert edge split if any
            insertEdgeSplit
            (
                edgeMidPoint,
                faceMidPointi,
                edgeMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointi,
                otherEdgeMidPointi,
                newFaceVerts
            );
           
            newFaceVerts.append(otherEdgeMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                otherEdgeMidPointi,  
                otherFaceMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(otherFaceMidPointi);
        }
        else if(is2D_ && (faceOrder != (mesh_.faceOwner()[facei] == celli)))
        {
            newFaceVerts.append(otherEdgeMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                otherEdgeMidPointi,
                edgeMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointi,
                faceMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(faceMidPointi);
            newFaceVerts.append(otherFaceMidPointi);
        }
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       

        // In 2D, i dont use anchorPointi but rather its projection on
        // the master patch, with is otherAnchorPointi.
        // We use it, together with otherAnchors, to set the own and nei of
        // new faces on 2D refinement

        face newFace;
        newFace.transfer(newFaceVerts);

        if(!is2D_) // PoliMi
        {
            label anchorCell0 = getAnchorCell
            (
                cellAnchorPoints,
                cellAddedCells,
                celli,
                facei,
                anchorPointi
            );
            label anchorCell1 = getAnchorCell
            (
                cellAnchorPoints,
                cellAddedCells,
                celli,
                facei,
                anchors.otherVertex(anchorPointi)
            );


            label own, nei;
            point ownPt, neiPt;

            if (anchorCell0 < anchorCell1)
            {
                own = anchorCell0;
                nei = anchorCell1;

                ownPt = mesh_.points()[anchorPointi];
                neiPt = mesh_.points()[anchors.otherVertex(anchorPointi)];

            }
            else
            {
                own = anchorCell1;
                nei = anchorCell0;
                newFace.flip();

                ownPt = mesh_.points()[anchors.otherVertex(anchorPointi)];
                neiPt = mesh_.points()[anchorPointi];
            }

            if (debug)
            {
                point ownPt, neiPt;

                if (anchorCell0 < anchorCell1)
                {
                    ownPt = mesh_.points()[anchorPointi];
                    neiPt = mesh_.points()[anchors.otherVertex(anchorPointi)];
                }
                else
                {
                    ownPt = mesh_.points()[anchors.otherVertex(anchorPointi)];
                    neiPt = mesh_.points()[anchorPointi];
                }

                checkInternalOrientation
                (
                    meshMod,
                    celli,
                    facei,
                    ownPt,
                    neiPt,
                    newFace
                );
            }

            return addInternalFace
            (
                meshMod,
                facei,
                anchorPointi,
                newFace,
                own,
                nei
            );
        }// end if is 2D
        else // is2D_
            // PoliMi
        {
            const edge& otherAnchors = midPointToAnchors[otherEdgeMidPointi];
 
            label anchorCell0 = getAnchorCell
            (
                cellAnchorPoints,
                cellAddedCells,
                celli,
                facei,
                otherAnchorPointi
            );
            label anchorCell1 = getAnchorCell
            (
                cellAnchorPoints,
                cellAddedCells,
                celli,
                facei,
                otherAnchors.otherVertex(otherAnchorPointi)
            );

            label own, nei;
            point ownPt, neiPt;

            if (anchorCell0 < anchorCell1)
            {
                own = anchorCell0;
                nei = anchorCell1;

                ownPt = mesh_.points()[otherAnchorPointi];
                neiPt = mesh_.points()[otherAnchors.otherVertex(otherAnchorPointi)];
            }
            else
            {
                own = anchorCell1;
                nei = anchorCell0;
                newFace.flip();

                ownPt = mesh_.points()[otherAnchors.otherVertex(otherAnchorPointi)];
                neiPt = mesh_.points()[otherAnchorPointi];
            }

            if (debug)
            {
                point ownPt, neiPt;

                if (anchorCell0 < anchorCell1)
                {
                    ownPt = mesh_.points()[otherAnchorPointi];
                    neiPt = mesh_.points()[otherAnchors.otherVertex(otherAnchorPointi)];
                }
                else
                {
                    ownPt = mesh_.points()[otherAnchors.otherVertex(otherAnchorPointi)];
                    neiPt = mesh_.points()[otherAnchorPointi];
                }

                checkInternalOrientation
                (
                    meshMod,
                    celli,
                    facei,
                    ownPt,
                    neiPt,
                    newFace
                );
            }

            return addInternalFace
            (
                meshMod,
                facei,
                otherAnchorPointi,
                newFace,
                own,
                nei
            );
        } // end !is2D_
    }
    else // not changed, 2 anchors or 2 faceMids
    {
        return -1;
    }
}


// Modified version for 2D
// Creates all the 4/12 internal faces for celli.
void Foam::hexRefPoliMi::createInternalFaces
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& faceAnchorLevel,
    const labelList& edgeMidPoint,
    const label celli,

    polyTopoChange& meshMod
) const
{
    // Find in every face the cellLevel+1 points (from edge subdivision)
    // and the anchor points.

    const cell& cFaces = mesh_.cells()[celli];
    const label cLevel = cellLevel_[celli];

    // From edge mid to anchor points
    Map<edge> midPointToAnchors(24);

    // From edge mid to face mids
    Map<edge> midPointToFaceMids(24);

    if(is2D_)// PoliMi
    {
        midPointToAnchors.resize(8);
        midPointToFaceMids.resize(8);
    }

    // Storage for on-the-fly addressing
    DynamicList<label> storage;


    // Running count of number of internal faces added so far.
    label nFacesAdded = 0;

    forAll(cFaces, i)
    {
        label facei = cFaces[i];

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // PoliMi
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Check if face is on an empty/wedge patch
       
        if(is2D_)
        {
            if(mesh_.isInternalFace(facei))
            {
                continue;
            }
            label patchi = mesh_.boundaryMesh().whichPatch(facei);
                    
            if (!(mesh_.boundaryMesh()[patchi].type() == "empty"
                || 
                mesh_.boundaryMesh()[patchi].type() == "wedge"))
            {
                // if not an empty/wedge face, we skip the face
                continue;
            }
        }
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        const face& f = mesh_.faces()[facei];
        const labelList& fEdges = mesh_.faceEdges(facei, storage);

        // We are on the celli side of face f. The face will have 1 or 4
        // cLevel points and lots of higher numbered ones.

        label faceMidPointi = -1;

        label nAnchors = countAnchors(f, cLevel);


        if (nAnchors == 1)
        {
            // Only one anchor point. So the other side of the face has already
            // been split using cLevel+1 and cLevel+2 points.

            // Find the one anchor.
            label anchorFp = -1;

            forAll(f, fp)
            {
                if (pointLevel_[f[fp]] <= cLevel)
                {
                    anchorFp = fp;
                    break;
                }
            }        
            // Now the face mid point is the second cLevel+1 point
            label edgeMid = findLevel
            (
                facei,
                f,
                f.fcIndex(anchorFp),
                true,
                cLevel+1
            );
            label faceMid = findLevel
            (
                facei,
                f,
                f.fcIndex(edgeMid),
                true,
                cLevel+1
            );

            faceMidPointi = f[faceMid];
        }
        else if (nAnchors == 4)
        {
            // There is no face middle yet but the face will be marked for
            // splitting.

            faceMidPointi = faceMidPoint[facei];
        }
        else
        {
            dumpCell(mesh_.faceOwner()[facei]);
            if (mesh_.isInternalFace(facei))
            {
                dumpCell(mesh_.faceNeighbour()[facei]);
            }

            FatalErrorInFunction
                << "nAnchors:" << nAnchors
                << " facei:" << facei
                << abort(FatalError);
        }



        // Now loop over all the anchors (might be just one) and store
        // the edge mids connected to it. storeMidPointInfo will collect
        // all the info and combine it all.
        forAll(f, fp0)
        {
            label point0 = f[fp0];

            if (pointLevel_[point0] <= cLevel)
            {
                // Anchor.

                // Walk forward
                // ~~~~~~~~~~~~
                // to cLevel+1 or edgeMidPoint of this level.

                label edgeMidPointi = -1;

                label fp1 = f.fcIndex(fp0);

                if (pointLevel_[f[fp1]] <= cLevel)
                {
                    // Anchor. Edge will be split.
                    label edgeI = fEdges[fp0];

                    edgeMidPointi = edgeMidPoint[edgeI];

                    if (edgeMidPointi == -1)
                    {
                        dumpCell(celli);

                        const labelList& cPoints = mesh_.cellPoints(celli);

                        FatalErrorInFunction
                            << "cell:" << celli << " cLevel:" << cLevel
                            << " cell points:" << cPoints
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, cPoints)()
                            << " face:" << facei
                            << " f:" << f
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, f)()
                            << " faceAnchorLevel:" << faceAnchorLevel[facei]
                            << " faceMidPoint:" << faceMidPoint[facei]
                            << " faceMidPointi:" << faceMidPointi
                            << " fp:" << fp0
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Search forward in face to clevel+1
                    label edgeMid = findLevel(facei, f, fp1, true, cLevel+1);

                    edgeMidPointi = f[edgeMid];
                }

                label newFacei = storeMidPointInfo
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    cellMidPoint,
                    edgeMidPoint,

                    celli,
                    facei,
                    true,                   // mid point after anchor
                    edgeMidPointi,          // edgemid
                    point0,                 // anchor
                    faceMidPointi,

                    midPointToAnchors,
                    midPointToFaceMids,
                    meshMod
                );

                if (newFacei != -1)
                {
                    nFacesAdded++;

                    if ((nFacesAdded == 12 && !is2D_)
                        ||
                        (nFacesAdded == 4 && is2D_))// PoliMi
                    {
                        break;
                    }
                }



                // Walk backward
                // ~~~~~~~~~~~~~

                label fpMin1 = f.rcIndex(fp0);

                if (pointLevel_[f[fpMin1]] <= cLevel)
                {
                    // Anchor. Edge will be split.
                    label edgeI = fEdges[fpMin1];

                    edgeMidPointi = edgeMidPoint[edgeI];

                    if (edgeMidPointi == -1)
                    {
                        dumpCell(celli);

                        const labelList& cPoints = mesh_.cellPoints(celli);

                        FatalErrorInFunction
                            << "cell:" << celli << " cLevel:" << cLevel
                            << " cell points:" << cPoints
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, cPoints)()
                            << " face:" << facei
                            << " f:" << f
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, f)()
                            << " faceAnchorLevel:" << faceAnchorLevel[facei]
                            << " faceMidPoint:" << faceMidPoint[facei]
                            << " faceMidPointi:" << faceMidPointi
                            << " fp:" << fp0
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Search back to clevel+1
                    label edgeMid = findLevel
                    (
                        facei,
                        f,
                        fpMin1,
                        false,
                        cLevel+1
                    );

                    edgeMidPointi = f[edgeMid];
                }

                newFacei = storeMidPointInfo
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    cellMidPoint,
                    edgeMidPoint,

                    celli,
                    facei,
                    false,                  // mid point before anchor
                    edgeMidPointi,          // edgemid
                    point0,                 // anchor
                    faceMidPointi,

                    midPointToAnchors,
                    midPointToFaceMids,
                    meshMod
                );

                if (newFacei != -1)
                {
                    nFacesAdded++;

                    if ((nFacesAdded == 12 && !is2D_)
                        ||
                        (nFacesAdded == 4 && is2D_)) // PoliMi
                    {
                        break;
                    }
                }
            }   // done anchor
        }   // done face

        if ((nFacesAdded == 12 && !is2D_)
            ||
            (nFacesAdded == 4 && is2D_)) // PoliMi
        {
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, read refinement data

Foam::hexRefPoliMi::hexRefPoliMi
    (const polyMesh& mesh, const bool readHistory /* = true*/)
:
    hexRef8PoliMi(mesh, readHistory),
    is2D_(is2D())
{};

Foam::hexRefPoliMi::hexRefPoliMi
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const refinementHistory& history,
    const scalar level0Edgei/* = -great*/
)
:
    hexRef8PoliMi(mesh, cellLevel, pointLevel, history, level0Edgei),
    is2D_(is2D())
{};
        
Foam::hexRefPoliMi::hexRefPoliMi
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const scalar level0Edge/* = -great*/
)
:
    hexRef8PoliMi(mesh, cellLevel, pointLevel, level0Edge),
    is2D_(is2D())
{};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// MODIFIED VERSION FOR 2D
// PoliMi

Foam::labelList Foam::hexRefPoliMi::getSplitPoints() const
{

    if (debug)
    {
        checkRefinementLevels(-1, labelList(0));
    }

    if (debug)
    {
        Pout<< "hexRefPoliMi::getSplitPoints :"
            << " Calculating unrefineable points" << endl;
    }


    if (!history_.active())
    {
        FatalErrorInFunction
            << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    // Master cell
    // -1 undetermined
    // -2 certainly not split point
    // >= label of master cell
    labelList splitMaster(mesh_.nPoints(), -1);
    labelList splitMasterLevel(mesh_.nPoints(), 0);

    // Unmark all with not 8 cells
    // const labelListList& pointCells = mesh_.pointCells();

    for (label pointi = 0; pointi < mesh_.nPoints(); pointi++)
    {
        const labelList& pCells = mesh_.pointCells(pointi);

        const label pCellsSize = pCells.size();
        if
        (
            (pCellsSize != 8 && !is2D_)
            ||
            (pCellsSize != 4 && is2D_) // PoliMi
        )
        {
            splitMaster[pointi] = -2;
        }
    }

    // Unmark all with different master cells
    const labelList& visibleCells = history_.visibleCells();
    //"Per cell in the current mesh (i.e. visible) either -1 (unrefined)
    // or an index into splitCells"

    forAll(visibleCells, celli)
    {
        const labelList& cPoints = mesh_.cellPoints(celli);

        if (visibleCells[celli] != -1 && history_.parentIndex(celli) >= 0)
        {
            label parentIndex = history_.parentIndex(celli);

            // Check same master.
            forAll(cPoints, i)
            {
                label pointi = cPoints[i];

                label masterCelli = splitMaster[pointi];

                if (masterCelli == -1)
                {
                    // First time visit of point. Store parent cell and
                    // level of the parent cell (with respect to celli). This
                    // is additional guarantee that we're referring to the
                    // same master at the same refinement level.

                    splitMaster[pointi] = parentIndex;
                    splitMasterLevel[pointi] = cellLevel_[celli] - 1;
                }
                else if (masterCelli == -2)
                {
                    // Already decided that point is not splitPoint
                }
                else if
                (
                    (masterCelli != parentIndex)
                 || (splitMasterLevel[pointi] != cellLevel_[celli] - 1)
                )
                {
                    // Different masters so point is on two refinement
                    // patterns
                    splitMaster[pointi] = -2;
                }
            }
        }
        else
        {
            // Either not visible or is unrefined cell
            forAll(cPoints, i)
            {
                label pointi = cPoints[i];

                splitMaster[pointi] = -2;
            }
        }
    }

    if(!is2D_)
    {
        // Unmark boundary faces
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            const face& f = mesh_.faces()[facei];

            forAll(f, fp)
            {
                splitMaster[f[fp]] = -2;
            }
        }
    }
    else// leave only points on master empty/wedge patch
    {
        // find master patch: the first appearing empty/wegde patch
        label masterPatch = -1;
        label slavePatch = -1;
        bool firstDone(false);

        for(int patchi = 0; patchi < mesh_.boundaryMesh().size(); patchi ++)
        {
            if (!firstDone)
            {
                if (mesh_.boundaryMesh()[patchi].type()=="empty" ||
                    mesh_.boundaryMesh()[patchi].type()=="wedge" )
                {
                    masterPatch = patchi;
                    firstDone = true;
                }
            }
            else
            {
                if (mesh_.boundaryMesh()[patchi].type()=="empty" ||
                    mesh_.boundaryMesh()[patchi].type()=="wedge" )
                {
                    slavePatch = patchi;
                    break;
                }
            }
        }

        if(masterPatch == -1 || slavePatch ==-1)
        {
            // abort, error
            FatalErrorInFunction
                << "Could not find empty/wedge patch in a 2D case"
                << abort(FatalError);
        }

        // delete points on salve patch
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            label patchFacei = mesh_.boundaryMesh().whichPatch(facei);

            if(patchFacei == slavePatch)
            {
                const face& f = mesh_.faces()[facei];

                forAll(f, fp)
                {
                    splitMaster[f[fp]] = -2;
                }
            }
        }
    }

    // Collect into labelList
    label nSplitPoints = 0;

    forAll(splitMaster, pointi)
    {
        if (splitMaster[pointi] >= 0)
        {
            nSplitPoints++;
        }
    }

    labelList splitPoints(nSplitPoints);
    nSplitPoints = 0;

    forAll(splitMaster, pointi)
    {
        if (splitMaster[pointi] >= 0)
        {
            splitPoints[nSplitPoints++] = pointi;
        }
    }

    return splitPoints;
}

// Top level driver to insert topo changes to do all refinement.
Foam::labelListList Foam::hexRefPoliMi::setRefinement
(
    const labelList& cellLabels,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();
        // Cannot call checkRefinementlevels since hanging points might
        // get triggered by the mesher after subsetting.
        // checkRefinementLevels(-1, labelList(0));
    }

    // Clear any saved point/cell data.
    savedPointLevel_.clear();
    savedCellLevel_.clear();


    // New point/cell level. Copy of pointLevel for existing points.
    DynamicList<label> newCellLevel(cellLevel_.size());
    forAll(cellLevel_, celli)
    {
        newCellLevel.append(cellLevel_[celli]);
    }
    DynamicList<label> newPointLevel(pointLevel_.size());
    forAll(pointLevel_, pointi)
    {
        newPointLevel.append(pointLevel_[pointi]);
    }


    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Allocating " << cellLabels.size() << " cell midpoints."
            << endl;
    }


    // Mid point per refined cell.
    // -1 : not refined
    // >=0: label of mid point.
    labelList cellMidPoint(mesh_.nCells(), -1);

    forAll(cellLabels, i)
    {
        label celli = cellLabels[i];

        if(!is2D_)
        {
            label anchorPointi = mesh_.faces()[mesh_.cells()[celli][0]][0];

            cellMidPoint[celli] = meshMod.setAction
            (
                polyAddPoint
                (
                    mesh_.cellCentres()[celli],     // point
                    anchorPointi,                   // master point
                    -1,                             // zone for point
                    true                            // supports a cell
                )
            );

            newPointLevel(cellMidPoint[celli]) = cellLevel_[celli]+1;
        }
        else // PoliMi
        {
            cellMidPoint[celli] = 1234567890;// only mark them as >0
        }
    }


    if (debug)
    {
        cellSet splitCells(mesh_, "splitCells", cellLabels.size());

        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                splitCells.insert(celli);
            }
        }

        Pout<< "hexRef8::setRefinement : Dumping " << splitCells.size()
            << " cells to split to cellSet " << splitCells.objectPath()
            << endl;

        splitCells.write();
    }

    // Split edges
    // ~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Allocating edge midpoints."
            << endl;
    }

    // Unrefined edges are ones between cellLevel or lower points.
    // If any cell using this edge gets split then the edge needs to be split.

    // -1  : no need to split edge
    // >=0 : label of introduced mid point
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PoliMi
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    boolList isRefinableFace(mesh_.nFaces(), false); // only for 2D
    // true for empty/wedge faces of cells marked for refinement
    //  mark all edges of these faces for refinement

    if(is2D_)
    {
        for (label facei = mesh_.nInternalFaces(); facei<mesh_.nFaces(); facei++)
            // TODO loop only over one empty/wedge patches?
        {
            const label& patchID = mesh_.boundaryMesh().whichPatch(facei);

            if((mesh_.boundaryMesh()[patchID].type()=="wedge")
                ||
                (mesh_.boundaryMesh()[patchID].type()=="empty"))
            {
                label own = mesh_.faceOwner()[facei];

                if(cellMidPoint[own] > 0)
                {
                    isRefinableFace[facei] = true;
                    const labelList& fEdges = mesh_.faceEdges(facei);

                    forAll(fEdges, k)
                    {
                        label edgei = fEdges[k];
                        const edge& e = mesh_.edges()[edgei];

                        if
                        (
                            pointLevel_[e[0]] <= cellLevel_[own]
                         && pointLevel_[e[1]] <= cellLevel_[own]
                        )
                        {
                            edgeMidPoint[edgei] = 1234567890;// mark for splitting
                        }
                    }
                }
            }
        }
    }
    else //if is 3D: original code
    {
        // Note: Loop over cells to be refined or edges?
        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                const labelList& cEdges = mesh_.cellEdges(celli);

                forAll(cEdges, i)
                {
                    label edgeI = cEdges[i];

                    const edge& e = mesh_.edges()[edgeI];

                    if
                    (
                        pointLevel_[e[0]] <= cellLevel_[celli]
                     && pointLevel_[e[1]] <= cellLevel_[celli]
                    )
                    {
                        edgeMidPoint[edgeI] = 12345;    // mark need for splitting
                    }
                }
            }
        }
    }

    // Synchronize edgeMidPoint across coupled patches. Take max so that
    // any split takes precedence.
    syncTools::syncEdgeList
    (
        mesh_,
        edgeMidPoint,
        maxEqOp<label>(),
        labelMin
    );


    // Introduce edge points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: calculate midpoints and sync.
        // This needs doing for if people do not write binary and we slowly
        // get differences.

        pointField edgeMids(mesh_.nEdges(), point(-great, -great, -great));

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split.
                edgeMids[edgeI] = mesh_.edges()[edgeI].centre(mesh_.points());
            }
        }
        syncTools::syncEdgePositions
        (
            mesh_,
            edgeMids,
            maxEqOp<vector>(),
            point(-great, -great, -great)
        );


        // Phase 2: introduce points at the synced locations.
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split. Replace edgeMidPoint with actual
                // point label.

                const edge& e = mesh_.edges()[edgeI];

                edgeMidPoint[edgeI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        edgeMids[edgeI],            // point
                        e[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                newPointLevel(edgeMidPoint[edgeI]) =
                    max
                    (
                        pointLevel_[e[0]],
                        pointLevel_[e[1]]
                    )
                  + 1;
            }
        }
    }

    if (debug)
    {
        OFstream str(mesh_.time().path()/"edgeMidPoint.obj");

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                const edge& e = mesh_.edges()[edgeI];

                meshTools::writeOBJ(str, e.centre(mesh_.points()));
            }
        }

        Pout<< "hexRef8::setRefinement :"
            << " Dumping edge centres to split to file " << str.name() << endl;
    }


    // Calculate face level
    // ~~~~~~~~~~~~~~~~~~~~
    // (after splitting)

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Allocating face midpoints."
            << endl;
    }

    // Face anchor level. There are guaranteed 4 points with level
    // <= anchorLevel. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());

    for (label facei = 0; facei < mesh_.nFaces(); facei++)
    {
        faceAnchorLevel[facei] = faceLevel(facei);
    }

    // -1  : no need to split face
    // >=0 : label of introduced mid point
    labelList faceMidPoint(mesh_.nFaces(), -1);


    // Internal faces: look at cells on both sides. Uniquely determined since
    // face itself guaranteed to be same level as most refined neighbour.
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceAnchorLevel[facei] >= 0)
        {
            label own = mesh_.faceOwner()[facei];
            label ownLevel = cellLevel_[own];
            label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

            label nei = mesh_.faceNeighbour()[facei];
            label neiLevel = cellLevel_[nei];
            label newNeiLevel = neiLevel + (cellMidPoint[nei] >= 0 ? 1 : 0);

            if
            (
                newOwnLevel > faceAnchorLevel[facei]
             || newNeiLevel > faceAnchorLevel[facei]
            )
            {
                faceMidPoint[facei] = 12345;    // mark to be split
            }
        }
    }

    // Coupled patches handled like internal faces except now all information
    // from neighbour comes from across processor.
    // Boundary faces are more complicated since the boundary face can
    // be more refined than its owner (or neighbour for coupled patches)
    // (does not happen if refining/unrefining only, but does e.g. when
    //  refinining and subsetting)

    {
        labelList newNeiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(newNeiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            label ownLevel = cellLevel_[own];
            label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

            newNeiLevel[i] = newOwnLevel;
        }

        // Swap.
        syncTools::swapBoundaryFaceList(mesh_, newNeiLevel);

        // So now we have information on the neighbour.

        forAll(newNeiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();

            if (faceAnchorLevel[facei] >= 0)
            {
                label own = mesh_.faceOwner()[facei];
                label ownLevel = cellLevel_[own];
                label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

                if
                (
                    newOwnLevel > faceAnchorLevel[facei]
                 || newNeiLevel[i] > faceAnchorLevel[facei]
                )
                {
                    faceMidPoint[facei] = 12345;    // mark to be split
                }
            }
        }
    }


    // Synchronize faceMidPoint across coupled patches. (logical or)
    syncTools::syncFaceList
    (
        mesh_,
        faceMidPoint,
        maxEqOp<label>()
    );



    // Introduce face points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: determine mid points and sync. See comment for edgeMids
        // above
        pointField bFaceMids
        (
            mesh_.nFaces()-mesh_.nInternalFaces(),
            point(-great, -great, -great)
        );

        forAll(bFaceMids, i)
        {
            label facei = i+mesh_.nInternalFaces();

            if (faceMidPoint[facei] >= 0)
            {
                bFaceMids[i] = mesh_.faceCentres()[facei];
            }
        }
        syncTools::syncBoundaryFacePositions
        (
            mesh_,
            bFaceMids,
            maxEqOp<vector>()
        );

        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0)
            {
                // Face marked to be split. Replace faceMidPoint with actual
                // point label.

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // PoliMi
                if(is2D_ && !isRefinableFace[facei])
                {
                    // skip internal faces in 2D
                    // They are marked as 1234567890 but the face mid point
                    // is not created
                    continue;
                }
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                const face& f = mesh_.faces()[facei];

                faceMidPoint[facei] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        (
                            facei < mesh_.nInternalFaces()
                          ? mesh_.faceCentres()[facei]
                          : bFaceMids[facei-mesh_.nInternalFaces()]
                        ),                          // point
                        f[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                // Determine the level of the corner points and midpoint will
                // be one higher.
                newPointLevel(faceMidPoint[facei]) = faceAnchorLevel[facei]+1;
            }
        }
    }

    if (debug)
    {
        faceSet splitFaces(mesh_, "splitFaces", cellLabels.size());

        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0)
            {
                splitFaces.insert(facei);
            }
        }

        Pout<< "hexRef8::setRefinement : Dumping " << splitFaces.size()
            << " faces to split to faceSet " << splitFaces.objectPath() << endl;

        splitFaces.write();
    }


    // Information complete
    // ~~~~~~~~~~~~~~~~~~~~
    // At this point we have all the information we need. We should no
    // longer reference the cellLabels to refine. All the information is:
    // - cellMidPoint >= 0 : cell needs to be split
    // - faceMidPoint >= 0 : face needs to be split
    // - edgeMidPoint >= 0 : edge needs to be split

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PoliMi
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // We would have the same for cellMidPoint and edgeMidPoint
    // For faceMidPoint we would have:
    //  - The label of center point for empty/wedge faces
    //  - 1234567890 for internal faces that need to be split in 2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Get the corner/anchor points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Finding cell anchorPoints (8 per cell)"
            << endl;
    }

    // There will always be 8 points on the hex that have were introduced
    // with the hex and will have the same or lower refinement level.

    // Per cell the 8 corner points.
    labelListList cellAnchorPoints(mesh_.nCells());

    {
        labelList nAnchorPoints(mesh_.nCells(), 0);

        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                cellAnchorPoints[celli].setSize(8);
            }
        }

        forAll(pointLevel_, pointi)
        {
            const labelList& pCells = mesh_.pointCells(pointi);

            forAll(pCells, pCelli)
            {
                label celli = pCells[pCelli];

                if
                (
                    cellMidPoint[celli] >= 0
                 && pointLevel_[pointi] <= cellLevel_[celli]
                )
                {
                    if (nAnchorPoints[celli] == 8)
                    {
                        dumpCell(celli);
                        FatalErrorInFunction
                            << "cell " << celli
                            << " of level " << cellLevel_[celli]
                            << " uses more than 8 points of equal or"
                            << " lower level" << nl
                            << "Points so far:" << cellAnchorPoints[celli]
                            << abort(FatalError);
                    }

                    cellAnchorPoints[celli][nAnchorPoints[celli]++] = pointi;
                }
            }
        }

        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                if (nAnchorPoints[celli] != 8)
                {
                    dumpCell(celli);

                    const labelList& cPoints = mesh_.cellPoints(celli);

                    FatalErrorInFunction
                        << "cell " << celli
                        << " of level " << cellLevel_[celli]
                        << " does not seem to have 8 points of equal or"
                        << " lower level" << endl
                        << "cellPoints:" << cPoints << endl
                        << "pointLevels:"
                        << UIndirectList<label>(pointLevel_, cPoints)() << endl
                        << abort(FatalError);
                }
            }
        }
    }


    // Add the cells
    // ~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Adding cells (1 per anchorPoint)"
            << endl;
    }

    // Per cell the 7 added cells (+ original cell)
    // PoliMi in 2D, we would add 3
    labelListList cellAddedCells(mesh_.nCells());

    forAll(cellAnchorPoints, celli)
    {
        const labelList& cAnchors = cellAnchorPoints[celli];

        if (cAnchors.size() == 8)
        {
            label nAddedCellsPlusOne = 8; // PoliMi
            if(is2D_)
            {
                nAddedCellsPlusOne = 4;
            }

            labelList& cAdded = cellAddedCells[celli];
            cAdded.setSize(nAddedCellsPlusOne);// PoliMi

            // Original cell at 0
            cAdded[0] = celli;

            // Update cell level
            newCellLevel[celli] = cellLevel_[celli]+1;


            for (label i = 1; i < nAddedCellsPlusOne; i++) // PoliMi
            {
                cAdded[i] = meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,                                 // master point
                        -1,                                 // master edge
                        -1,                                 // master face
                        celli,                              // master cell
                        mesh_.cellZones().whichZone(celli)  // zone for cell
                    )
                );

                newCellLevel(cAdded[i]) = cellLevel_[celli]+1;
            }
        }
    }


    // Faces
    // ~~~~~
    // 1. existing faces that get split (into four always)
    // 2. existing faces that do not get split but only edges get split
    // 3. existing faces that do not get split but get new owner/neighbour
    // 4. new internal faces inside split cells.

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Marking faces to be handled"
            << endl;
    }

    // Get all affected faces.
    PackedBoolList affectedFace(mesh_.nFaces());

    {
        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                const cell& cFaces = mesh_.cells()[celli];

                forAll(cFaces, i)
                {
                    affectedFace.set(cFaces[i]);
                }
            }
        }

        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0)
            {
                affectedFace.set(facei);
            }
        }

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                const labelList& eFaces = mesh_.edgeFaces(edgeI);

                forAll(eFaces, i)
                {
                    affectedFace.set(eFaces[i]);
                }
            }
        }
    }


    // 1. Faces that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef8::setRefinement : Splitting faces" << endl;
    }

    forAll(faceMidPoint, facei)
    {
        if (faceMidPoint[facei] >= 0 && affectedFace.get(facei))
        {
            // Face needs to be split and hasn't yet been done in some way
            // (affectedFace - is impossible since this is first change but
            //  just for completeness)

            const face& f = mesh_.faces()[facei];

            // Has original facei been used (three faces added, original gets
            // modified)
            bool modifiedFace = false;

            label anchorLevel = faceAnchorLevel[facei];

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // PoliMi
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if(is2D_ && !isRefinableFace[facei])// face to split in 2
            {
                // Maybe in a completely separated function??
                // face to split in 2
                face newFace(2);
                
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    label nextPointi = f[f.fcIndex(fp)];
                    label edgei = meshTools::findEdge(mesh_, pointi, nextPointi);

                    if(edgeMidPoint[edgei] >= 0)
                    {
                        DynamicList<label> faceVerts(4);
                        faceVerts.append(pointi);

                        label pointj = f[f.rcIndex(fp)];


                        walkFaceToMid
                        (
                            edgeMidPoint,
                            anchorLevel,
                            facei,
                            fp,
                            faceVerts
                        );

                        walkFaceFromMid
                        (
                            edgeMidPoint,
                            anchorLevel,
                            facei,
                            f.rcIndex(fp),
                            faceVerts
                        );

                        faceVerts.append(pointj);
                        
                        // Convert dynamiclist to face.
                        newFace.transfer(faceVerts);
                        
                        // Get new owner/neighbour
                        label own, nei;
                        getFaceNeighbours
                        (
                            cellAnchorPoints,
                            cellAddedCells,
                            facei,
                            pointi,          // Anchor point

                            own,
                            nei
                        );

                        if (debug)
                        {
                            if (mesh_.isInternalFace(facei))
                            {
                                label oldOwn = mesh_.faceOwner()[facei];
                                label oldNei = mesh_.faceNeighbour()[facei];

                                checkInternalOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.cellCentres()[oldNei],
                                    newFace
                                );
                            }
                            else
                            {
                                label oldOwn = mesh_.faceOwner()[facei];

                                checkBoundaryOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.faceCentres()[facei],
                                    newFace
                                );
                            }
                        }


                        if (!modifiedFace)
                        {
                            modifiedFace = true;

                            modFace(meshMod, facei, newFace, own, nei);
                        }
                        else
                        {
                            addFace(meshMod, facei, newFace, own, nei);
                        }
                    }
                }
            }
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            else
            {
                // original code to split face in 4
                face newFace(4);

                forAll(f, fp)
                {
                    label pointi = f[fp];

                    if (pointLevel_[pointi] <= anchorLevel)
                    {
                        // point is anchor. Start collecting face.

                        DynamicList<label> faceVerts(4);

                        faceVerts.append(pointi);

                        // Walk forward to mid point.
                        // - if next is +2 midpoint is +1
                        // - if next is +1 it is midpoint
                        // - if next is +0 there has to be edgeMidPoint

                        walkFaceToMid
                        (
                            edgeMidPoint,
                            anchorLevel,
                            facei,
                            fp,
                            faceVerts
                        );

                        faceVerts.append(faceMidPoint[facei]);

                        walkFaceFromMid
                        (
                            edgeMidPoint,
                            anchorLevel,
                            facei,
                            fp,
                            faceVerts
                        );

                        // Convert dynamiclist to face.
                        newFace.transfer(faceVerts);

                        // Pout<< "Split face:" << facei << " verts:" << f
                        //    << " into quad:" << newFace << endl;

                        // Get new owner/neighbour
                        label own, nei;
                        getFaceNeighbours
                        (
                            cellAnchorPoints,
                            cellAddedCells,
                            facei,
                            pointi,          // Anchor point

                            own,
                            nei
                        );


                        if (debug)
                        {
                            if (mesh_.isInternalFace(facei))
                            {
                                label oldOwn = mesh_.faceOwner()[facei];
                                label oldNei = mesh_.faceNeighbour()[facei];

                                checkInternalOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.cellCentres()[oldNei],
                                    newFace
                                );
                            }
                            else
                            {
                                label oldOwn = mesh_.faceOwner()[facei];

                                checkBoundaryOrientation
                                (
                                    meshMod,
                                    oldOwn,
                                    facei,
                                    mesh_.cellCentres()[oldOwn],
                                    mesh_.faceCentres()[facei],
                                    newFace
                                );
                            }
                        }


                        if (!modifiedFace)
                        {
                            modifiedFace = true;

                            modFace(meshMod, facei, newFace, own, nei);
                        }
                        else
                        {
                            addFace(meshMod, facei, newFace, own, nei);
                        }
                    }
                }
            } // end if is not face to split in 2

            // Mark face as having been handled
            affectedFace.unset(facei);
        }
    }

    // 2. faces that do not get split but use edges that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Adding edge splits to unsplit faces"
            << endl;
    }

    DynamicList<label> eFacesStorage;
    DynamicList<label> fEdgesStorage;

    forAll(edgeMidPoint, edgeI)
    {
        if (edgeMidPoint[edgeI] >= 0)
        {
            // Split edge. Check that face not already handled above.

            const labelList& eFaces = mesh_.edgeFaces(edgeI, eFacesStorage);

            forAll(eFaces, i)
            {
                label facei = eFaces[i];

                if (faceMidPoint[facei] < 0 && affectedFace.get(facei))
                {
                    // Unsplit face. Add edge splits to face.

                    const face& f = mesh_.faces()[facei];
                    const labelList& fEdges = mesh_.faceEdges
                    (
                        facei,
                        fEdgesStorage
                    );

                    DynamicList<label> newFaceVerts(f.size());

                    forAll(f, fp)
                    {
                        newFaceVerts.append(f[fp]);

                        label edgeI = fEdges[fp];

                        if (edgeMidPoint[edgeI] >= 0)
                        {
                            newFaceVerts.append(edgeMidPoint[edgeI]);
                        }
                    }

                    face newFace;
                    newFace.transfer(newFaceVerts);

                    // The point with the lowest level should be an anchor
                    // point of the neighbouring cells.
                    label anchorFp = findMinLevel(f);

                    label own, nei;
                    getFaceNeighbours
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        facei,
                        f[anchorFp],          // Anchor point

                        own,
                        nei
                    );


                    if (debug)
                    {
                        if (mesh_.isInternalFace(facei))
                        {
                            label oldOwn = mesh_.faceOwner()[facei];
                            label oldNei = mesh_.faceNeighbour()[facei];

                            checkInternalOrientation
                            (
                                meshMod,
                                oldOwn,
                                facei,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.cellCentres()[oldNei],
                                newFace
                            );
                        }
                        else
                        {
                            label oldOwn = mesh_.faceOwner()[facei];

                            checkBoundaryOrientation
                            (
                                meshMod,
                                oldOwn,
                                facei,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.faceCentres()[facei],
                                newFace
                            );
                        }
                    }

                    modFace(meshMod, facei, newFace, own, nei);

                    // Mark face as having been handled
                    affectedFace.unset(facei);
                }
            }
        }
    }


    // 3. faces that do not get split but whose owner/neighbour change
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Changing owner/neighbour for otherwise unaffected faces"
            << endl;
    }

    forAll(affectedFace, facei)
    {
        if (affectedFace.get(facei))
        {
            const face& f = mesh_.faces()[facei];

            // The point with the lowest level should be an anchor
            // point of the neighbouring cells.
            label anchorFp = findMinLevel(f);

            label own, nei;
            getFaceNeighbours
            (
                cellAnchorPoints,
                cellAddedCells,
                facei,
                f[anchorFp],          // Anchor point

                own,
                nei
            );

            modFace(meshMod, facei, f, own, nei);

            // Mark face as having been handled
            affectedFace.unset(facei);
        }
    }


    // 4. new internal faces inside split cells.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    // This is the hard one. We have to find the splitting points between
    // the anchor points. But the edges between the anchor points might have
    // been split (into two,three or four edges).

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Create new internal faces for split cells"
            << endl;
    }

    forAll(cellMidPoint, celli)
    {
        if (cellMidPoint[celli] >= 0)
        {
            createInternalFaces
            (
                cellAnchorPoints,
                cellAddedCells,
                cellMidPoint,
                faceMidPoint,
                faceAnchorLevel,
                edgeMidPoint,
                celli,
                meshMod
            );
        }
    }

    // Extend pointLevels and cellLevels for the new cells. Could also be done
    // in updateMesh but saves passing cellAddedCells out of this routine.

    // Check
    if (debug)
    {
        label minPointi = labelMax;
        label maxPointi = labelMin;

        forAll(cellMidPoint, celli)
        {
            if (cellMidPoint[celli] >= 0)
            {
                minPointi = min(minPointi, cellMidPoint[celli]);
                maxPointi = max(maxPointi, cellMidPoint[celli]);
            }
        }
        forAll(faceMidPoint, facei)
        {
            if (faceMidPoint[facei] >= 0)
            {
                minPointi = min(minPointi, faceMidPoint[facei]);
                maxPointi = max(maxPointi, faceMidPoint[facei]);
            }
        }
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                minPointi = min(minPointi, edgeMidPoint[edgeI]);
                maxPointi = max(maxPointi, edgeMidPoint[edgeI]);
            }
        }

        if (minPointi != labelMax && minPointi != mesh_.nPoints())
        {
            FatalErrorInFunction
                << "Added point labels not consecutive to existing mesh points."
                << nl
                << "mesh_.nPoints():" << mesh_.nPoints()
                << " minPointi:" << minPointi
                << " maxPointi:" << maxPointi
                << abort(FatalError);
        }
    }

    pointLevel_.transfer(newPointLevel);
    cellLevel_.transfer(newCellLevel);

    // Mark files as changed
    setInstance(mesh_.facesInstance());


    // Update the live split cells tree.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // New unrefinement structure
    if (history_.active())
    {
        if (debug)
        {
            Pout<< "hexRef8::setRefinement :"
                << " Updating refinement history to " << cellLevel_.size()
                << " cells" << endl;
        }

        // Extend refinement history for new cells
        history_.resize(cellLevel_.size());

        forAll(cellAddedCells, celli)
        {
            const labelList& addedCells = cellAddedCells[celli];

            if (addedCells.size())
            {
                // Cell was split.
                history_.storeSplit(celli, addedCells);
            }
        }
    }

    // Compact cellAddedCells.

    labelListList refinedCells(cellLabels.size());

    forAll(cellLabels, i)
    {
        label celli = cellLabels[i];

        refinedCells[i].transfer(cellAddedCells[celli]);
    }

    return refinedCells;
}

void Foam::hexRefPoliMi::setUnrefinement
(
    const labelList& splitPointLabels,
    polyTopoChange& meshMod
)
{

    if (!history_.active())
    {
        FatalErrorInFunction
            << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "hexRefPoliMi::setUnrefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();

        forAll(cellLevel_, celli)
        {
            if (cellLevel_[celli] < 0)
            {
                FatalErrorInFunction
                    << "Illegal cell level " << cellLevel_[celli]
                    << " for cell " << celli
                    << abort(FatalError);
            }
        }


        // Write to sets.
        pointSet pSet(mesh_, "splitPoints", splitPointLabels);
        pSet.write();

        cellSet cSet(mesh_, "splitPointCells", splitPointLabels.size());

        forAll(splitPointLabels, i)
        {
            const labelList& pCells = mesh_.pointCells(splitPointLabels[i]);

            forAll(pCells, j)
            {
                cSet.insert(pCells[j]);
            }
        }
        cSet.write();

        Pout<< "hexRefPoliMi::setRefinement : Dumping " << pSet.size()
            << " points and "
            << cSet.size() << " cells for unrefinement to" << nl
            << "    pointSet " << pSet.objectPath() << nl
            << "    cellSet " << cSet.objectPath()
            << endl;
    }


    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    {
        labelHashSet splitFaces(12*splitPointLabels.size());
        if(!is2D_) // PoliMi
        {
            forAll(splitPointLabels, i)
            {
                const labelList& pFaces = mesh_.pointFaces()[splitPointLabels[i]];

                forAll(pFaces, j)
                {
                    splitFaces.insert(pFaces[j]);
                }
            }
        }
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // PoliMi
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        else
        {
            splitFaces.resize(4*splitPointLabels.size());

            forAll(splitPointLabels, i)
            {
                const labelList& pFaces = mesh_.pointFaces()[splitPointLabels[i]];

                forAll(pFaces, j)
                {
                    if(mesh_.isInternalFace(pFaces[j]))
                        // only remove internal faces
                    {
                        splitFaces.insert(pFaces[j]);
                    }
                }
            }
        }
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Check with faceRemover what faces will get removed. Note that this
        // can be more (but never less) than splitFaces provided.
        faceRemover_.compatibleRemoves
        (
            splitFaces.toc(),   // pierced faces
            cellRegion,         // per cell -1 or region it is merged into
            cellRegionMaster,   // per region the master cell
            facesToRemove       // new faces to be removed.
        );

        if (facesToRemove.size() != splitFaces.size())
        {
            FatalErrorInFunction
                << "Initial set of split points to unrefine does not"
                << " seem to be consistent or not mid points of refined cells"
                << abort(FatalError);
        }
    }

    // Redo the region master so it is consistent with our master.
    // This will guarantee that the new cell (for which faceRemover uses
    // the region master) is already compatible with our refinement structure.

    forAll(splitPointLabels, i)
    {
        label pointi = splitPointLabels[i];

        // Get original cell label

        const labelList& pCells = mesh_.pointCells(pointi);

        // Check
        if ((!is2D_ && pCells.size() != 8) ||
                (is2D_ && pCells.size() != 4)) // PoliMi
        {
            FatalErrorInFunction
                << "splitPoint " << pointi
                << " should have 8 cells in 3D or 4 in 2D using it. It has " << pCells
                << abort(FatalError);
        }


        // Check that the lowest numbered pCells is the master of the region
        // (should be guaranteed by directRemoveFaces)
        // if (debug)
        {
            label masterCelli = min(pCells);

            forAll(pCells, j)
            {
                label celli = pCells[j];

                label region = cellRegion[celli];

                if (region == -1)
                {
                    FatalErrorInFunction
                        << "Initial set of split points to unrefine does not"
                        << " seem to be consistent or not mid points"
                        << " of refined cells" << nl
                        << "cell:" << celli << " on splitPoint " << pointi
                        << " has no region to be merged into"
                        << abort(FatalError);
                }

                if (masterCelli != cellRegionMaster[region])
                {
                    FatalErrorInFunction
                        << "cell:" << celli << " on splitPoint:" << pointi
                        << " in region " << region
                        << " has master:" << cellRegionMaster[region]
                        << " which is not the lowest numbered cell"
                        << " among the pointCells:" << pCells
                        << abort(FatalError);
                }
            }
        }
    }

    // Insert all commands to combine cells. Never fails so don't have to
    // test for success.
    faceRemover_.setRefinement
    (
        facesToRemove,
        cellRegion,
        cellRegionMaster,
        meshMod
    );
    // Remove the 8 cells that originated from merging around the split point
    // and adapt cell levels (not that pointLevels stay the same since points
    // either get removed or stay at the same position.
    forAll(splitPointLabels, i)
    {
        label pointi = splitPointLabels[i];

        const labelList& pCells = mesh_.pointCells(pointi);

        label masterCelli = min(pCells);

        forAll(pCells, j)
        {
            cellLevel_[pCells[j]]--;
        }

        history_.combineCells(masterCelli, pCells);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // history_.updateMesh will take care of truncating.
}


// ************************************************************************* //
