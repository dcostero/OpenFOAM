/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::hexRefPoliMi

Description
    Refinement of (split) hexes using polyTopoChange.

SourceFiles
    hexRefPoliMi.C

\*---------------------------------------------------------------------------*/

#ifndef hexRefPoliMi_H
#define hexRefPoliMi_H

#include "labelIOList.H"
#include "face.H"
#include "HashSet.H"
#include "DynamicList.H"
#include "primitivePatch.H"
#include "removeFaces.H"
//#include "refinementHistory.H"
#include "PackedBoolList.H"
#include "uniformDimensionedFields.H"
#include "cellShapeList.H"
#include "hexRef8.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyPatch;
class polyTopoChange;
class mapPolyMesh;
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                           Class hexRefPoliMi Declaration
\*---------------------------------------------------------------------------*/

class hexRefPoliMi : public hexRef8
{
    // Private Data

        const bool is2D_;

    // Private Member Functions
       
        bool is2D();

        //- Get cell added to point of celli (if any)
        label getAnchorCell
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const label celli,
            const label facei,
            const label pointi
        ) const;

        //- Store in maps correspondence from midpoint to anchors and faces.
        label storeMidPointInfo
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& edgeMidPoint,
            const label celli,
            const label facei,
            const bool faceOrder,
            const label midPointi,
            const label anchorPointi,
            const label faceMidPointi,

            Map<edge>& midPointToAnchors,
            Map<edge>& midPointToFaceMids,
            polyTopoChange& meshMod
        ) const;

        //- Create all internal faces to split celli into 8.
        void createInternalFaces
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& faceMidPoint,
            const labelList& faceAnchorLevel,
            const labelList& edgeMidPoint,
            const label celli,
            polyTopoChange& meshMod
        ) const;

        // Cellshape recognition

public:

    //- Runtime type information
    ClassName("hexRefPoliMi");


    // Constructors

        //- Construct from mesh, read_if_present refinement data
        //  (from write below). If readHistory is true does read_if_present
        //  of refinement history. If false clears all history
        hexRefPoliMi(const polyMesh& mesh, const bool readHistory = true);

        //- Construct from mesh and un/refinement data and optional size of
        //  starting cells
        hexRefPoliMi
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const refinementHistory& history,
            const scalar level0Edge = -great
        );

        //- Construct from mesh and refinement data and optional size of
        //  starting cells
        hexRefPoliMi
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const scalar level0Edge = -great
        );


  

    // Member Functions

        // Access

       // Refinement

            //- Insert refinement. All selected cells will be split into 8.
            //  Returns per element in cells the 8 cells they were split into.
            //  Guarantees that the 0th element is the original cell label.
            //  Mapping:
            //  -split cells: 7 new ones get added from original
            //  -split faces: original gets modified; 3 new ones get added
            //               from original
            //  -added internal faces: added from original cell face(if
            //   that was internal) or created out-of-nothing (so will not
            //   get mapped!). Note: could make this inflate from point but
            //   that will allocate interpolation.
            //  -points added to split edge: added from edge start()
            //  -midpoints added: added from cellPoints[0].
            labelListList setRefinement
            (
                const labelList& cells,
                polyTopoChange&
            );

        // Unrefinement (undoing refinement, not arbitrary coarsening)

            virtual labelList getSplitPoints() const;

            //- Remove some refinement. Needs to be supplied output of
            //  consistentUnrefinement. Only call if undoable set.
            //  All 8 pointCells of a split point will be combined into
            //  the lowest numbered cell of those 8.
            virtual void setUnrefinement
            (
                const labelList& splitPointLabels,
                polyTopoChange&
            );


            virtual void writeObj(){};

    
    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const hexRefPoliMi&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
