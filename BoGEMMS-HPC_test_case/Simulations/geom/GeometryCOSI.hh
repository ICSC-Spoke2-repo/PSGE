/***************************************************************************
                          GeometryCOSI.hh  -  description
                             -------------------
 begin                : 2022
 Authors              : V. Fioretti (INAF OAS Bologna) valentina.fioretti@inaf.it

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software for non commercial purpose              *
 *   and for public research institutes; you can redistribute it and/or    *
 *   modify it under the terms of the GNU General Public License.          *
 *   For commercial purpose see appropriate license terms                  *
 *                                                                         *
 ***************************************************************************/
// ********************************************************************
//
// Copyright 2025 Valentina Fioretti
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// **********************************************************************

#ifndef GeometryCOSI_H
#define GeometryCOSI_H 1
#include "globals.hh"


class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4Region;
class G4UserLimits;

#include "G4VisAttributes.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Material.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Element.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "MaterialsDefinition.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//Factory
#include "GeoClassFactory.hh"


class GeometryCOSI: public G4VUserDetectorConstruction {
public:

    	GeometryCOSI();
    	~GeometryCOSI();

    	G4VPhysicalVolume* Construct();
        void SetPhysicalWorld(G4VPhysicalVolume* World_phys);
    	void DefineMaterials();
    	G4VPhysicalVolume* ConstructGeometry(G4VPhysicalVolume* World_phys);

protected:

    
    	void DefineSensitiveDetector();
        MaterialsDefinition* materials;

    	// World
    	//Physical Volumes
    	G4VPhysicalVolume* World_phys;
    	G4LogicalVolume* World_log;
    
        //Materials
    
        G4Material* chamber_mat;
        G4Material* bgo_mat;
        G4Material* csi_mat;        
        G4Material* refl1_mat;
        G4Material* refl2_mat;
        G4Material* glue_mat;
        G4Material* pmt_mat;
        G4Material* SiPad_mat;
        G4Material* boro_glass_mat;
        G4Material* bialkali_mat;
        G4Material* air_mat;
private:
    
    G4UserLimits* target_limit;            // pointer to user step limits
    static DerivedRegister<GeometryCOSI> reg;
};

#endif

