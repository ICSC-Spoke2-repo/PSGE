/***************************************************************************
                          GeometryCOSI.cc  -  description
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

#include "GeometryCOSI.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Polyhedra.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "THELGlobalMemory.hh"
#include "globals.hh"
#include "G4NistManager.hh"
#include "G4PVReplica.hh"
#include "MaterialsDefinition.hh"
#include "G4PVParameterised.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UserLimits.hh"

// Regions
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

#include <vector>

//GDML
//#include "G4Writer/G4GDMLWriter.h"

// CADMesh
#include "CADMesh.hh"


GeometryCOSI::GeometryCOSI()
: G4VUserDetectorConstruction(), target_limit(NULL)
{

	World_phys = 0;
	World_log = 0;
	materials = new MaterialsDefinition;

}

GeometryCOSI::~GeometryCOSI() {

    delete target_limit;

}

DerivedRegister<GeometryCOSI> GeometryCOSI::reg("GeometryCOSI");

G4VPhysicalVolume* GeometryCOSI::Construct() {

    //DefineMaterials();
    World_phys = gm.ConstructWorld();
    World_log = gm.World_log;
    World_phys = ConstructGeometry(World_phys);
    return World_phys;
	
}


G4VPhysicalVolume* GeometryCOSI::ConstructGeometry(G4VPhysicalVolume* World_phys) {
	

    // --------------------- PARAMETERS -----------------------

    // CLAIRE CsI experiment      
    //%%%%%%%%% LABORATORY CHAMBER
    
    // vacuum or not vacuum in the laboratory
    G4int vacuum_flag = 0;
    gm.config->readInto(vacuum_flag, "GEOM.COSI.VACUUM");
    G4cout << "GEOM.COSI.VACUUM: " << vacuum_flag << G4endl;
        
    if (vacuum_flag){
        chamber_mat = materials->GetMaterial(12);
    }
    if (!vacuum_flag){
        chamber_mat = materials->GetMaterial(45);
    }
        
        
    //chamber
    G4double chamber_side = 2.0*m;
    G4Box* chamber = new G4Box("chamber",
                            chamber_side/2.,
                            chamber_side/2.,
                            chamber_side/2.);
        
    G4LogicalVolume* log_chamber = new G4LogicalVolume(chamber, chamber_mat, "log_chamber");
    G4VPhysicalVolume* phys_chamber = new G4PVPlacement(0, G4ThreeVector(0.0,0.0,0.0), log_chamber, "phys_chamber", World_log, false, 0);
        
    log_chamber->SetVisAttributes (&G4VisAttributes::GetInvisible);
    //gm.AddXYZDetector( log_chamber);
        

        
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// CsI set-up parameters


    G4int CsI_type = 0; //CsI type: 1 = CLAIRE experiment
    gm.config->readInto(CsI_type, "GEOM.COSI.CsI.TYPE");
    G4cout << "GEOM.COSI.CsI.TYPE: " << CsI_type << G4endl;
        
    G4int ReflLayer1_mat = 0; // reflective layer 1 mat 
    gm.config->readInto(ReflLayer1_mat, "GEOM.COSI.2.REFL.LAYER1.MAT");
    G4cout << "GEOM.COSI.2.REFL.LAYER1.MAT: " << ReflLayer1_mat << G4endl;

    G4int ReflLayer2_mat = 0; // reflective layer 2 mat 
    gm.config->readInto(ReflLayer2_mat, "GEOM.COSI.2.REFL.LAYER2.MAT");
    G4cout << "GEOM.COSI.2.REFL.LAYER2.MAT: " << ReflLayer2_mat << G4endl;
        
    G4int setReflSurfaceType = 1; // 1: ground (2 deg.), 2:
    gm.config->readInto(setReflSurfaceType, "PHYS.COSI.2.OPTSURFACE.WRAPPER");
    G4cout << "PHYS.COSI.2.OPTSURFACE.WRAPPER: " << setReflSurfaceType << G4endl;
        
    // %%%%%%%%%%%%%%%%%%%%%%%%%% Materials
    csi_mat = materials->GetMaterial(14);
    refl1_mat = materials->GetMaterial(ReflLayer1_mat);
    refl2_mat = materials->GetMaterial(ReflLayer2_mat);
    glue_mat = materials->GetMaterial(25);
    pmt_mat = materials->GetMaterial(12); 
    SiPad_mat = materials->GetMaterial(63);
    boro_glass_mat = materials->GetMaterial(65);
    bialkali_mat = materials->GetMaterial(66);
    air_mat = materials->GetMaterial(11);

        
    // %%%%%%%%%%%%%%%%%%%%%%%%%% Material optical properties
    
    // chamber mat
    G4MaterialPropertiesTable* chamber_MPT = new G4MaterialPropertiesTable();
        
    std::vector<G4double> chamber_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
    std::vector<G4double> chamber_RIND    = { 1., 1., 1., 1., 1.};
    chamber_MPT->AddProperty("RINDEX", chamber_Energy, chamber_RIND);
    chamber_MPT->DumpTable();
        
    chamber_mat->SetMaterialPropertiesTable(chamber_MPT);
    
    // CsI
    std::vector<G4double> CsI_energy_emission = {1.76013558*eV, 1.7756688*eV,  1.79147862*eV, 1.80989526*eV, 1.8168995*eV,  1.82159919*eV,
        1.83345547*eV, 1.8430522*eV,  1.85274993*eV, 1.86501654*eV, 1.87494741*eV, 1.88246524*eV,
        1.8900436*eV,  1.90024352*eV, 1.91055412*eV, 1.92360078*eV, 1.93682686*eV, 1.94485017*eV,
        1.95023608*eV, 1.96109788*eV, 1.96657427*eV, 1.98040002*eV, 1.99442154*eV, 2.01438858*eV,
        2.02889728*eV, 2.03770326*eV, 2.04956417*eV, 2.05855085*eV, 2.07370509*eV, 2.07370509*eV,
        2.0829052*eV,  2.08908409*eV, 2.10469291*eV, 2.11734891*eV, 2.13338457*eV, 2.14312309*eV,
        2.14638905*eV, 2.15955304*eV, 2.16953249*eV, 2.1796046*eV,  2.1966009*eV,  2.21386434*eV,
        2.24205746*eV, 2.26732211*eV, 2.30065419*eV, 2.33111634*eV, 2.35844021*eV, 2.37434338*eV,
        2.38237567*eV, 2.39452651*eV, 2.41092176*eV, 2.42754308*eV, 2.44439517*eV, 2.44439517*eV,
        2.46579218*eV, 2.47445624*eV, 2.4875671*eV,  2.50526591*eV, 2.51421009*eV, 2.51421009*eV,
        2.52774676*eV, 2.54602403*eV, 2.55990639*eV, 2.56456754*eV, 2.56456754*eV, 2.56456754*eV,
        2.5692457*eV,  2.58338315*eV, 2.60247689*eV, 2.61698345*eV, 2.62185498*eV, 2.63165264*eV,
        2.64648721*eV, 2.66652877*eV, 2.68176027*eV, 2.70753651*eV, 2.72851704*eV, 2.76060465*eV,
        2.79345596*eV, 2.83278458*eV, 2.867387*eV,   2.90284521*eV, 2.93919135*eV, 2.9701824*eV,
        3.01468428*eV, 3.06053997*eV, 3.10781221*eV, 3.15656767*eV, 3.21419556*eV, 3.27396675*eV,
        3.32027462*eV, 3.39224589*eV, 3.46740644*eV, 3.52820765*eV, 3.59117923*eV, 3.6564395*eV};

    std::vector<G4double> CsI_emission = {0.15189873, 0.16455696, 0.17721519, 0.19240506, 0.21265823, 0.23544304,
        0.2556962,  0.27341772, 0.28860759, 0.30886076, 0.33164557, 0.34936709,
        0.3721519,  0.39240506, 0.41518987, 0.43544304, 0.4556962,  0.47594937,
        0.49367089, 0.5164557,  0.53670886, 0.55949367, 0.57974684, 0.59493671,
        0.61012658, 0.63291139, 0.6556962,  0.67341772, 0.69367089, 0.71392405,
        0.73924051, 0.75949367, 0.77974684, 0.8,        0.82278481, 0.84303797,
        0.86582278, 0.88607595, 0.90379747, 0.92658228, 0.94683544, 0.96708861,
        0.97721519, 0.98227848, 0.97974684, 0.9721519,  0.95949367, 0.93670886,
        0.91392405, 0.89367089, 0.87341772, 0.85063291, 0.80759494, 0.78481013,
        0.76455696, 0.7443038,  0.72151899, 0.70126582, 0.67848101, 0.65316456,
        0.63291139, 0.61265823, 0.59240506, 0.54683544, 0.52405063, 0.50126582,
        0.56708861, 0.48101266, 0.46075949, 0.43797468, 0.41518987, 0.39240506,
        0.3721519,  0.34936709, 0.32911392, 0.30886076, 0.28860759, 0.27088608,
        0.2556962,  0.23797468, 0.22278481, 0.2,        0.18481013, 0.16962025,
        0.15189873, 0.13670886, 0.12405063, 0.10886076, 0.09367089, 0.07848101,
        0.06329114, 0.05822785, 0.04303797, 0.03037975, 0.01772152, 0.00759494};
    
    std::vector<G4double> CsI_energy_rindex = {0.0185051*eV,  0.01956821*eV, 0.02069508*eV, 0.021886*eV,   0.02314433*eV, 0.02447379*eV,
        0.02587856*eV, 0.02736958*eV, 0.02894122*eV, 0.03060583*eV, 0.0323634*eV,  0.03423087*eV,
        0.03619977*eV, 0.03827854*eV, 0.04047803*eV, 0.04281222*eV, 0.04526623*eV, 0.04787035*eV,
        0.05062646*eV, 0.05353376*eV, 0.05661379*eV, 0.05986683*eV, 0.06332186*eV, 0.06694611*eV,
        0.07080765*eV, 0.07486969*eV, 0.07917254*eV, 0.08371654*eV, 0.08856014*eV, 0.09364365*eV,
        0.09902891*eV, 0.10471638*eV, 0.11070018*eV, 0.11707667*eV, 0.12386034*eV, 0.13096461*eV,
        0.13848341*eV, 0.14644956*eV, 0.15488345*eV, 0.16378362*eV, 0.17321067*eV, 0.18316472*eV,
        0.19369505*eV, 0.20483099*eV, 0.21660412*eV, 0.22904895*eV, 0.24225127*eV, 0.2561657*eV,
        0.27088529*eV, 0.28646996*eV, 0.30291766*eV, 0.32037261*eV, 0.33875464*eV, 0.3582323*eV,
        0.37880904*eV, 0.40059515*eV, 0.42373274*eV, 0.44808167*eV, 0.47376461*eV, 0.50114874*eV,
        0.529847*eV,   0.56025395*eV, 0.59265869*eV, 0.62649923*eV, 0.66266274*eV, 0.70087167*eV,
        0.74108905*eV, 0.78371807*eV, 0.82877138*eV, 0.87621342*eV, 0.92663825*eV, 0.98011224*eV,
        1.03665718*eV, 1.09623518*eV, 1.15873083*eV, 1.2251403*eV,  1.2962279*eV,  1.37074846*eV,
        1.44959895*eV, 1.53294014*eV, 1.6211323*eV,  1.71438327*eV, 1.81289952*eV, 1.9171826*eV,
        2.02721057*eV, 2.14394256*eV, 2.26703599*eV, 2.3976832*eV,  2.53546418*eV, 2.68131917*eV,
        2.83522064*eV, 2.99840867*eV, 3.17095137*eV, 3.35273657*eV, 3.54645877*eV, 3.75027823*eV,
        3.96495678*eV, 4.19290492*eV, 4.43434186*eV, 4.6892662*eV,  4.95936794*eV,};
    
    
    std::vector<G4double> CsI_rindex = {1.52755217, 1.55614332, 1.5803263,  1.60087763, 1.61845679, 1.63357901,
        1.64665327, 1.65805045, 1.6679514,  1.67662173, 1.68421077, 1.69090889,
        1.69678608, 1.70196171, 1.70653687, 1.71059758, 1.71417313, 1.71735514,
        1.72018116, 1.72268528, 1.72491567, 1.72689754, 1.72866922, 1.73023492,
        1.73164128, 1.73288895, 1.73400457, 1.73499964, 1.73589557, 1.73669036,
        1.73740295, 1.7380404,  1.73860933, 1.73912414, 1.73958934, 1.740004,
        1.74037857, 1.74071806, 1.74102637, 1.74130649, 1.74156306, 1.74179852,
        1.74201635, 1.74221922, 1.74240963, 1.74258997, 1.74276313, 1.74293015,
        1.74309381, 1.74325628, 1.74341909, 1.74358521, 1.74375545, 1.743933,
        1.74411958, 1.74431796, 1.74453136, 1.74476053, 1.74500877, 1.7452821,
        1.74557926, 1.74590714, 1.74627225, 1.74667171, 1.74711993, 1.74761822,
        1.74817089, 1.74878897, 1.74947904, 1.75024735, 1.75111163, 1.75208269,
        1.75317157, 1.75438901, 1.75574477, 1.75727534, 1.75901842, 1.76096436,
        1.76315877, 1.7656334,  1.76843079, 1.77159491, 1.77517618, 1.77924469,
        1.78386064, 1.78913892, 1.79515375, 1.80207378, 1.81001264, 1.81919193,
        1.82982176, 1.84226344, 1.85688665, 1.87415272, 1.89499245, 1.9201518,
        1.95104475, 1.99012362, 2.04094337, 2.1096908,  2.20937662};
    
    
    G4MaterialPropertiesTable* CsI_MPT = new G4MaterialPropertiesTable();
    //https://cds.cern.ch/record/2268719/files/1-s2.0-S0168900215009754-main.pdf
    
    // property independent of energy
    CsI_MPT->AddConstProperty("SCINTILLATIONYIELD", 54000.0/MeV);
    
    // decay time
    CsI_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 1000.*ns);
    // resolution scale
    CsI_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
    // rise time
    //CsI_MPT->AddConstProperty("SCINTILLATIONRISETIME1", 0.*ns);
    
    // properties that depend on energy
    CsI_MPT->AddProperty("SCINTILLATIONCOMPONENT1", CsI_energy_emission, CsI_emission);
    
    // properties that depend on energy
    CsI_MPT->AddProperty("RINDEX", CsI_energy_rindex, CsI_rindex, false, true);
    
    G4double CsI_absl_length = 0.; // absorption length
    gm.config->readInto(CsI_absl_length, "PHYS.COSI.2.CsI.ABSL");
    G4cout << "PHYS.COSI.2.CsI.ABSL: " << CsI_absl_length << G4endl;

    std::vector<G4double> CsI_Energy_abs  = { 0.001 * eV, 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
    std::vector<G4double> CsI_abs_length  = { CsI_absl_length*mm, CsI_absl_length*mm, CsI_absl_length*mm, CsI_absl_length*mm, CsI_absl_length*mm, CsI_absl_length*mm };
    CsI_MPT->AddProperty("ABSLENGTH", CsI_Energy_abs, CsI_abs_length, false, true);

    CsI_MPT->DumpTable();
    csi_mat->SetMaterialPropertiesTable(CsI_MPT);
    

    // reflLayer1(teflon)
    
    std::vector<G4double> reflLayer1_Energy  = { 0.001 * eV, 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
    std::vector<G4double> reflLayer1_RIND    = { 1.38, 1.38, 1.38, 1.38, 1.38, 1.38};

    G4MaterialPropertiesTable* reflLayer1_MPT = new G4MaterialPropertiesTable();
    reflLayer1_MPT->AddProperty("RINDEX", reflLayer1_Energy, reflLayer1_RIND, false, true);
    
    if (setReflSurfaceType == 7) {
        std::vector<G4double> spec_spike_Tef = {0., 0., 0., 0., 0., 0.};        // Direct reflection

        std::vector<G4double> spec_lobe_Tef = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2};    // Direct spread reflection -> lobe shape

        std::vector<G4double> back_scattering_Tef = {0., 0., 0., 0., 0., 0.};        // Back reflection

        reflLayer1_MPT->AddProperty("SPECULARSPIKECONSTANT", reflLayer1_Energy, spec_spike_Tef);
        reflLayer1_MPT->AddProperty("SPECULARLOBECONSTANT", reflLayer1_Energy, spec_lobe_Tef);
        reflLayer1_MPT->AddProperty("BACKSCATTERCONSTANT", reflLayer1_Energy, back_scattering_Tef);
    }

    reflLayer1_MPT->DumpTable();
    refl1_mat->SetMaterialPropertiesTable(reflLayer1_MPT);
    
    // reflLayer2(Aluminum)
    
    std::vector<G4double> reflLayer2_Energy  = {1.24000318*eV,1.25000452*eV,1.26000202*eV,1.27000459*eV,1.27999544*eV,1.28999707
        *eV,1.30000627*eV,1.3100059*eV,1.32000595*eV,1.33000288*eV,1.33999307*eV,1.35000216
        *eV,1.35999779*eV,1.37000628*eV,1.37999419*eV,1.39000413*eV,1.40000224*eV,1.41000089
        *eV,1.419997*eV,1.4300039*eV,1.4400023*eV,1.45000583*eV,1.45999456*eV,1.46999986
        *eV,1.48000189*eV,1.4899977*eV,1.5000024*eV,1.50999523*eV,1.51999164*eV,1.530008
        *eV,1.5400042*eV,1.54999623*eV,1.56000099*eV,1.56999656*eV,1.58000023*eV,1.59000985
        *eV,1.60000256*eV,1.60999621*eV,1.62000965*eV,1.6299984*eV,1.64000262*eV,1.64999865
        *eV,1.66000614*eV,1.67000079*eV,1.68000269*eV,1.69000993*eV,1.69999724*eV,1.71000894
        *eV,1.71999609*eV,1.73000402*eV,1.740007*eV,1.7500028*eV,1.75998919*eV,1.76998913
        *eV,1.78000113*eV,1.78999781*eV,1.80000288*eV,1.8099883*eV,1.82000497*eV,1.82999806
        *eV,1.83999226*eV,1.85001341*eV,1.86000478*eV,1.86999183*eV,1.88000119*eV,1.89000302
        *eV,1.89999538*eV,1.91000568*eV,1.92000307*eV,1.9299855*eV,1.94001155*eV,1.94998897
        *eV,1.96000756*eV,1.97000442*eV,1.98000892*eV,1.98998778*eV,2.0000032*eV,2.0099896
        *eV,2.0200104*eV,2.02999866*eV,2.0399855*eV,2.05000328*eV,2.05998303*eV,2.06999129
        *eV,2.07999259*eV,2.08998531*eV,2.10000336*eV,2.11001018*eV,2.12000408*eV,2.12998331
        *eV,2.13998306*eV,2.15000257*eV,2.16000346*eV,2.16998387*eV,2.18001861*eV,2.18999185
        *eV,2.20001772*eV,2.21001762*eV,2.21998959*eV,2.23001184*eV,2.24000359*eV,2.2500036
        *eV,2.26001091*eV,2.26998294*eV,2.28000144*eV,2.28998187*eV,2.30000739*eV,2.30999196
        *eV,2.32002018*eV,2.33000448*eV,2.33998676*eV,2.3500104*eV,2.3599855*eV,2.37000035
        *eV,2.38000918*eV,2.39001076*eV,2.40000384*eV,2.40998714*eV,2.42000661*eV,2.43001447
        *eV,2.44000942*eV,2.44999009*eV,2.46000394*eV,2.47000156*eV,2.47998157*eV,2.48999254
        *eV,2.49998384*eV,2.51000483*eV,2.52000403*eV,2.52997997*eV,2.53998317*eV,2.55001334
        *eV,2.56001731*eV,2.56999354*eV,2.57999414*eV,2.59001877*eV,2.60001255*eV,2.60997386
        *eV,2.6200118*eV,2.63001566*eV,2.63998378*eV,2.65002775*eV,2.65997722*eV,2.67000169
        *eV,2.67998613*eV,2.68998717*eV,2.70000432*eV,2.70997789*eV,2.72002541*eV,2.73002749
        *eV,2.73998229*eV,2.75000995*eV,2.75998839*eV,2.76997762*eV,2.77997709*eV,2.78998624
        *eV,2.80000448*eV,2.81003124*eV,2.82000178*eV,2.82997874*eV,2.84002654*eV,2.85001491
        *eV,2.86000781*eV,2.87000459*eV,2.88000461*eV,2.89000719*eV,2.90001166*eV,2.91001733
        *eV,2.92002351*eV,2.9300295*eV,2.94003458*eV,2.94996784*eV,2.95996845*eV,2.96996595
        *eV,2.98003121*eV,2.9900207*eV,3.0000048*eV,3.00998273*eV,3.02002724*eV,3.02999092
        *eV,3.04002056*eV,3.04996675*eV,3.05997824*eV,3.06997966*eV,3.07997015*eV,3.09002588
        *eV,3.09999246*eV,3.11002354*eV,3.11996272*eV,3.12996563*eV,3.14003288*eV,3.15000504
        *eV,3.1599602*eV,3.16997848*eV,3.17997893*eV,3.18996059*eV,3.20000512*eV,3.21002999
        *eV,3.22003424*eV,3.23001689*eV,3.23997696*eV,3.24999865*eV,3.2599968*eV,3.26997042
        *eV,3.28000525*eV,3.29001455*eV,3.2999973*eV,3.31004081*eV,3.31996783*eV,3.33004401
        *eV,3.34000157*eV,3.35001887*eV,3.36000538*eV,3.36996*eV,3.37997379*eV,3.38995457
        *eV,3.39999447*eV,3.41000023*eV,3.41997072*eV,3.42999968*eV,3.43999219*eV,3.45004309
        *eV,3.45995977*eV,3.47003074*eV,3.47996515*eV,3.48995661*eV,3.5000056*eV,3.51001326
        *eV,3.51997838*eV,3.53000024*eV,3.53997826*eV,3.55001284*eV,3.56000225*eV,3.57004804
        *eV,3.58004731*eV,3.5899988*eV,3.60000576*eV,3.60996356*eV,3.6199766*eV,3.63004534
        *eV,3.6399565*eV,3.65002939*eV,3.66005014*eV,3.67001742*eV,3.68003913*eV,3.69000591
        *eV,3.70002681*eV,3.70999128*eV,3.72000955*eV,3.72996987*eV,3.73998366*eV,3.75005137
        *eV,3.75994537*eV,3.77000634*eV,3.78000605*eV,3.7899431*eV,3.80004899*eV,3.80997475
        *eV,3.8199525*eV,3.82998265*eV,3.83994668*eV,3.84996269*eV,3.86003108*eV,3.87003148
        *eV,3.8799624*eV,3.88994442*eV,3.89997793*eV,3.90994003*eV,3.91995316*eV,3.9300177
        *eV,3.94000885*eV,3.95005092*eV,3.96001784*eV,3.97003517*eV,3.97997555*eV,3.98996584
        *eV,4.0000064*eV,4.00996793*eV,4.0199792*eV,4.03004058*eV,4.0400208*eV,4.05005058
        *eV,4.05999733*eV,4.06999305*eV,4.08003812*eV,4.08999797*eV,4.10000656*eV,4.11006426
        *eV,4.12003451*eV,4.13005325*eV,4.13998258*eV,4.14995978*eV,4.15998518*eV,4.17005914
        *eV,4.18004108*eV,4.18992932*eV,4.20000672*eV,4.20998976*eV,4.22002037*eV,4.22995457
        *eV,4.23993566*eV,4.24996395*eV,4.2600398*eV,4.27001648*eV,4.28003999*eV,4.28996223
        *eV,4.29993058*eV,4.30994537*eV,4.32000691*eV,4.32996432*eV,4.33996774*eV,4.35001749
        *eV,4.35996056*eV,4.36994919*eV,4.37998369*eV,4.39006439*eV,4.40003543*eV,4.41005188
        *eV,4.41995645*eV,4.4300639*eV,4.44005867*eV,4.44993893*eV,4.46002369*eV,4.46999309
        *eV,4.48000717*eV,4.49006622*eV,4.5000072*eV,4.5099923*eV,4.52002182*eV,4.52993052
        *eV,4.54004901*eV,4.55004582*eV,4.55991903*eV,4.57000363*eV,4.57996374*eV,4.58996736
        *eV,4.60001478*eV,4.60993487*eV,4.62007*eV,4.63007687*eV,4.63995354*eV,4.65004682
        *eV,4.66000896*eV,4.67001388*eV,4.68006185*eV,4.68997573*eV,4.69993171*eV,4.70993004
        *eV,4.71997101*eV,4.73005488*eV,4.7400007*eV,4.74998845*eV,4.76001837*eV,4.77009074
        *eV,4.78002153*eV,4.78999376*eV,4.80000768*eV,4.81006356*eV,4.81997428*eV,4.82992592
        *eV,4.83991874*eV,4.849953*eV,4.86002895*eV,4.86995555*eV,4.87992279*eV,4.88993092
        *eV,4.89998018*eV,4.91007083*eV,4.92000787*eV,4.92998523*eV,4.94000313*eV,4.95006182
        *eV,4.95996313*eV,4.96990413*eV,4.98008509*eV,4.98990616*eV,4.99996767*eV,5.01006984
        *eV,5.02000965*eV,5.02998898*eV,5.04000807*eV,5.05006714*eV,5.05995994*eV,5.0700989
        *eV,5.08007041*eV,5.09008122*eV,5.09992178*eV,5.11001106*eV,5.11992891*eV,5.13009758
        *eV,5.14009363*eV,5.14991478*eV,5.15998828*eV,5.17010126*eV,5.18003754*eV,5.19001207
        *eV,5.2000251*eV,5.21007683*eV,5.21994773*eV,5.23007671*eV,5.2400236*eV,5.2500084
        *eV,5.26003133*eV,5.2700926*eV,5.27996757*eV,5.29010532*eV,5.30005551*eV,5.31004319
        *eV,5.32006859*eV,5.32990278*eV,5.34000338*eV,5.34991148*eV,5.36008813*eV,5.37007096
        *eV,5.38009106*eV,5.38991429*eV,5.40000864*eV,5.40990481*eV,5.42007425*eV,5.43004417
        *eV,5.44005083*eV,5.45009444*eV,5.45993476*eV,5.47005199*eV,5.47996457*eV,5.48991314
        *eV,5.4998979*eV};
    
    std::vector<G4double> reflLayer2_RIND = { 9.45280e-01, 9.70930e-01, 9.96370e-01, 1.02147e+00, 1.04652e+00, 1.07150e+00
        , 1.09684e+00, 1.12312e+00, 1.15085e+00, 1.18042e+00, 1.21266e+00, 1.24832e+00
        , 1.28752e+00, 1.33076e+00, 1.37779e+00, 1.42868e+00, 1.48243e+00, 1.53824e+00
        , 1.59457e+00, 1.64985e+00, 1.70258e+00, 1.75087e+00, 1.79304e+00, 1.82839e+00
        , 1.85601e+00, 1.87539e+00, 1.88656e+00, 1.88979e+00, 1.88587e+00, 1.87531e+00
        , 1.85909e+00, 1.83807e+00, 1.81296e+00, 1.78488e+00, 1.75427e+00, 1.72192e+00
        , 1.68859e+00, 1.65413e+00, 1.61955e+00, 1.58489e+00, 1.55061e+00, 1.51655e+00
        , 1.48310e+00, 1.45052e+00, 1.41881e+00, 1.38753e+00, 1.35763e+00, 1.32829e+00
        , 1.30012e+00, 1.27267e+00, 1.24623e+00, 1.22072e+00, 1.19586e+00, 1.17210e+00
        , 1.14916e+00, 1.12692e+00, 1.10544e+00, 1.08469e+00, 1.06475e+00, 1.04542e+00
        , 1.02678e+00, 1.00885e+00, 9.91420e-01, 9.74560e-01, 9.58270e-01, 9.42440e-01
        , 9.27230e-01, 9.12500e-01, 8.98150e-01, 8.84370e-01, 8.70860e-01, 8.57830e-01
        , 8.45220e-01, 8.32900e-01, 8.20990e-01, 8.09400e-01, 7.98150e-01, 7.87190e-01
        , 7.76570e-01, 7.66180e-01, 7.56120e-01, 7.46350e-01, 7.36770e-01, 7.27510e-01
        , 7.18450e-01, 7.09570e-01, 7.00960e-01, 6.92560e-01, 6.84350e-01, 6.76330e-01
        , 6.68510e-01, 6.60860e-01, 6.53420e-01, 6.46100e-01, 6.38990e-01, 6.32010e-01
        , 6.25190e-01, 6.18510e-01, 6.11960e-01, 6.05560e-01, 5.99290e-01, 5.93190e-01
        , 5.87160e-01, 5.81280e-01, 5.75550e-01, 5.69880e-01, 5.64350e-01, 5.58890e-01
        , 5.53550e-01, 5.48320e-01, 5.43180e-01, 5.38160e-01, 5.33180e-01, 5.28340e-01
        , 5.23550e-01, 5.18860e-01, 5.14270e-01, 5.09760e-01, 5.05300e-01, 5.00910e-01
        , 4.96640e-01, 4.92410e-01, 4.88240e-01, 4.84170e-01, 4.80110e-01, 4.76180e-01
        , 4.72280e-01, 4.68430e-01, 4.64670e-01, 4.60960e-01, 4.57320e-01, 4.53700e-01
        , 4.50150e-01, 4.46640e-01, 4.43220e-01, 4.39820e-01, 4.36470e-01, 4.33180e-01
        , 4.29920e-01, 4.26720e-01, 4.23590e-01, 4.20470e-01, 4.17390e-01, 4.14370e-01
        , 4.11380e-01, 4.08440e-01, 4.05540e-01, 4.02660e-01, 3.99840e-01, 3.97040e-01
        , 3.94300e-01, 3.91580e-01, 3.88900e-01, 3.86260e-01, 3.83640e-01, 3.81050e-01
        , 3.78500e-01, 3.75970e-01, 3.73500e-01, 3.71040e-01, 3.68600e-01, 3.66200e-01
        , 3.63830e-01, 3.61490e-01, 3.59180e-01, 3.56890e-01, 3.54640e-01, 3.52390e-01
        , 3.50200e-01, 3.48020e-01, 3.45870e-01, 3.43720e-01, 3.41610e-01, 3.39540e-01
        , 3.37480e-01, 3.35430e-01, 3.33410e-01, 3.31420e-01, 3.29440e-01, 3.27500e-01
        , 3.25570e-01, 3.23650e-01, 3.21760e-01, 3.19910e-01, 3.18060e-01, 3.16220e-01
        , 3.14410e-01, 3.12620e-01, 3.10850e-01, 3.09100e-01, 3.07350e-01, 3.05640e-01
        , 3.03930e-01, 3.02250e-01, 3.00580e-01, 2.98930e-01, 2.97290e-01, 2.95680e-01
        , 2.94060e-01, 2.92470e-01, 2.90900e-01, 2.89350e-01, 2.87800e-01, 2.86280e-01
        , 2.84760e-01, 2.83260e-01, 2.81780e-01, 2.80310e-01, 2.78840e-01, 2.77390e-01
        , 2.75970e-01, 2.74550e-01, 2.73140e-01, 2.71750e-01, 2.70380e-01, 2.69010e-01
        , 2.67660e-01, 2.66320e-01, 2.64990e-01, 2.63660e-01, 2.62360e-01, 2.61070e-01
        , 2.59790e-01, 2.58500e-01, 2.57250e-01, 2.56000e-01, 2.54750e-01, 2.53510e-01
        , 2.52300e-01, 2.51090e-01, 2.49900e-01, 2.48710e-01, 2.47530e-01, 2.46360e-01
        , 2.45210e-01, 2.44050e-01, 2.42910e-01, 2.41780e-01, 2.40660e-01, 2.39560e-01
        , 2.38450e-01, 2.37350e-01, 2.36270e-01, 2.35200e-01, 2.34130e-01, 2.33060e-01
        , 2.32020e-01, 2.30970e-01, 2.29940e-01, 2.28930e-01, 2.27910e-01, 2.26890e-01
        , 2.25900e-01, 2.24900e-01, 2.23910e-01, 2.22930e-01, 2.21970e-01, 2.21000e-01
        , 2.20040e-01, 2.19080e-01, 2.18150e-01, 2.17220e-01, 2.16290e-01, 2.15370e-01
        , 2.14460e-01, 2.13550e-01, 2.12640e-01, 2.11750e-01, 2.10880e-01, 2.09990e-01
        , 2.09110e-01, 2.08250e-01, 2.07390e-01, 2.06530e-01, 2.05690e-01, 2.04850e-01
        , 2.04010e-01, 2.03170e-01, 2.02360e-01, 2.01540e-01, 2.00730e-01, 1.99910e-01
        , 1.99120e-01, 1.98320e-01, 1.97540e-01, 1.96750e-01, 1.95970e-01, 1.95200e-01
        , 1.94420e-01, 1.93680e-01, 1.92900e-01, 1.92150e-01, 1.91420e-01, 1.90660e-01
        , 1.89940e-01, 1.89210e-01, 1.88490e-01, 1.87770e-01, 1.87040e-01, 1.86340e-01
        , 1.85630e-01, 1.84920e-01, 1.84220e-01, 1.83530e-01, 1.82840e-01, 1.82160e-01
        , 1.81480e-01, 1.80810e-01, 1.80140e-01, 1.79470e-01, 1.78810e-01, 1.78160e-01
        , 1.77500e-01, 1.76860e-01, 1.76210e-01, 1.75570e-01, 1.74930e-01, 1.74300e-01
        , 1.73670e-01, 1.73050e-01, 1.72430e-01, 1.71820e-01, 1.71210e-01, 1.70600e-01
        , 1.70000e-01, 1.69390e-01, 1.68800e-01, 1.68210e-01, 1.67620e-01, 1.67030e-01
        , 1.66450e-01, 1.65880e-01, 1.65300e-01, 1.64730e-01, 1.64160e-01, 1.63600e-01
        , 1.63040e-01, 1.62480e-01, 1.61930e-01, 1.61380e-01, 1.60830e-01, 1.60290e-01
        , 1.59750e-01, 1.59220e-01, 1.58680e-01, 1.58160e-01, 1.57630e-01, 1.57100e-01
        , 1.56580e-01, 1.56060e-01, 1.55550e-01, 1.55040e-01, 1.54530e-01, 1.54030e-01
        , 1.53530e-01, 1.53030e-01, 1.52530e-01, 1.52040e-01, 1.51550e-01, 1.51060e-01
        , 1.50580e-01, 1.50100e-01, 1.49620e-01, 1.49140e-01, 1.48670e-01, 1.48200e-01
        , 1.47730e-01, 1.47260e-01, 1.46800e-01, 1.46340e-01, 1.45890e-01, 1.45440e-01
        , 1.44980e-01, 1.44530e-01, 1.44090e-01, 1.43640e-01, 1.43200e-01, 1.42760e-01
        , 1.42330e-01, 1.41890e-01, 1.41460e-01, 1.41030e-01, 1.40610e-01, 1.40180e-01
        , 1.39760e-01, 1.39340e-01, 1.38920e-01, 1.38510e-01, 1.38100e-01, 1.37690e-01
        , 1.37280e-01, 1.36870e-01, 1.36470e-01, 1.36070e-01, 1.35670e-01, 1.35270e-01
        , 1.34880e-01, 1.34480e-01, 1.34100e-01, 1.33710e-01, 1.33320e-01, 1.32940e-01
        , 1.32560e-01, 1.32180e-01, 1.31800e-01, 1.31430e-01, 1.31050e-01, 1.30680e-01
        , 1.30310e-01, 1.29950e-01, 1.29580e-01, 1.29220e-01, 1.28850e-01, 1.28490e-01
        , 1.28140e-01, 1.27780e-01, 1.27430e-01, 1.27080e-01, 1.26730e-01, 1.26380e-01
        , 1.26030e-01, 1.25690e-01, 1.25340e-01, 1.25000e-01, 1.24660e-01, 1.24310e-01
        , 1.23970e-01};
    
    std::vector<G4double> reflLayer2_KIND = {7.69295, 7.63287, 7.57296, 7.51368, 7.45468, 7.39616, 7.33792, 7.27925, 7.22077
        , 7.16271, 7.10504, 7.04828, 6.99369, 6.9419, 6.89437, 6.8518, 6.81583, 6.78728
        , 6.76703, 6.75544, 6.75239, 6.75748, 6.76979, 6.78793, 6.81056, 6.83616, 6.86302
        , 6.88965, 6.91538, 6.93857, 6.95904, 6.97607, 6.98941, 6.99877, 7.00432, 7.00614
        , 7.00438, 6.99915, 6.99085, 6.97968, 6.96595, 6.94967, 6.93142, 6.9113, 6.88956
        , 6.86617, 6.84172, 6.81589, 6.78935, 6.76178, 6.73351, 6.70464, 6.67511, 6.64536
        , 6.61525, 6.58475, 6.55393, 6.52308, 6.49213, 6.46103, 6.4298, 6.39889, 6.36765
        , 6.3367, 6.30565, 6.27465, 6.24386, 6.21318, 6.18269, 6.15241, 6.12212, 6.09212
        , 6.06229, 6.03255, 6.00308, 5.9737, 5.94467, 5.91571, 5.88711, 5.85859, 5.83024
        , 5.80242, 5.7745, 5.74694, 5.71954, 5.69232, 5.66536, 5.63867, 5.61207, 5.58582
        , 5.55976, 5.53378, 5.50817, 5.48265, 5.4575, 5.43235, 5.40757, 5.38289, 5.3584
        , 5.33411, 5.3101, 5.28629, 5.26268, 5.23917, 5.21606, 5.19295, 5.17015, 5.14736
        , 5.12479, 5.10251, 5.08036, 5.05842, 5.03659, 5.01489, 4.99351, 4.97214, 4.95111
        , 4.9302, 4.90941, 4.88875, 4.86833, 4.84814, 4.82787, 4.80795, 4.78806, 4.76841
        , 4.7489, 4.72953, 4.7103, 4.69121, 4.67238, 4.65348, 4.63494, 4.61633, 4.59798
        , 4.57979, 4.56163, 4.54364, 4.5258, 4.50811, 4.49059, 4.47312, 4.45569, 4.43866
        , 4.42156, 4.40452, 4.38776, 4.37094, 4.35441, 4.33793, 4.32163, 4.30538, 4.2892
        , 4.27319, 4.25725, 4.24149, 4.22579, 4.21015, 4.1947, 4.17931, 4.164, 4.14887
        , 4.13381, 4.11882, 4.1039, 4.08918, 4.07453, 4.05995, 4.04557, 4.03114, 4.01679
        , 4.00252, 3.98845, 3.97445, 3.96054, 3.94671, 3.93283, 3.91916, 3.90558, 3.89208
        , 3.87879, 3.86546, 3.85221, 3.83905, 3.82599, 3.813, 3.79998, 3.78718, 3.77447
        , 3.76185, 3.74919, 3.73663, 3.72416, 3.71179, 3.69951, 3.68732, 3.67524, 3.66312
        , 3.65109, 3.63917, 3.62721, 3.61535, 3.60372, 3.59193, 3.58038, 3.56879, 3.55744
        , 3.54592, 3.53465, 3.52334, 3.51214, 3.50105, 3.48993, 3.47891, 3.46787, 3.45693
        , 3.44611, 3.43539, 3.42465, 3.41402, 3.40336, 3.39281, 3.38238, 3.37192, 3.36157
        , 3.3512, 3.34094, 3.33081, 3.32064, 3.3106, 3.30052, 3.29042, 3.28044, 3.27058
        , 3.26069, 3.25093, 3.2413, 3.23163, 3.22194, 3.21238, 3.20279, 3.19333, 3.18384
        , 3.17448, 3.16511, 3.15581, 3.14656, 3.13736, 3.12822, 3.11914, 3.11005, 3.10101
        , 3.09206, 3.08311, 3.07419, 3.06532, 3.05655, 3.04778, 3.03905, 3.03036, 3.02172
        , 3.01312, 3.00456, 2.99606, 2.98757, 2.97916, 2.97076, 2.96242, 2.95411, 2.94584
        , 2.93762, 2.92943, 2.92127, 2.91317, 2.90508, 2.89704, 2.88902, 2.88106, 2.87311
        , 2.86521, 2.85732, 2.8495, 2.84173, 2.834, 2.82629, 2.81858, 2.81091, 2.80327
        , 2.79571, 2.78817, 2.78065, 2.77312, 2.76569, 2.75827, 2.75089, 2.74356, 2.73624
        , 2.72893, 2.72163, 2.71445, 2.70726, 2.70009, 2.69292, 2.68581, 2.67876, 2.67172
        , 2.66472, 2.65772, 2.65076, 2.64382, 2.63691, 2.63006, 2.62324, 2.61642, 2.60965
        , 2.60287, 2.59613, 2.5894, 2.58274, 2.57609, 2.56948, 2.56291, 2.55635, 2.5498
        , 2.54327, 2.53677, 2.53033, 2.5239, 2.5175, 2.5111, 2.50472, 2.4984, 2.49209
        , 2.48581, 2.47955, 2.47329, 2.46708, 2.46089, 2.45471, 2.44858, 2.44245, 2.43637
        , 2.43031, 2.42426, 2.41826, 2.41226, 2.40627, 2.40033, 2.39439, 2.38846, 2.38258
        , 2.37673, 2.37088, 2.36506, 2.35925, 2.35345, 2.34769, 2.34197, 2.33623, 2.33056
        , 2.32488, 2.31923, 2.31361, 2.30797, 2.30237, 2.29679, 2.29125, 2.28569, 2.28021
        , 2.27473, 2.26927, 2.26381, 2.25839, 2.25295, 2.24755, 2.24217, 2.23681, 2.23148
        , 2.22615, 2.22084, 2.21554, 2.21029, 2.20506, 2.19981, 2.19461, 2.18941, 2.18424
        , 2.17908, 2.17393, 2.16881, 2.16369, 2.15862, 2.15354, 2.14849, 2.14346, 2.13843
        , 2.13342, 2.12842, 2.12346, 2.11851, 2.11355, 2.10864, 2.10373, 2.09885, 2.09396
        , 2.08911, 2.08427, 2.07942, 2.0746, 2.06983, 2.06504, 2.06027, 2.05553, 2.05078
        , 2.04605, 2.04133, 2.03667, 2.03198, 2.02733, 2.02269, 2.01806, 2.01343, 2.00883
        , 2.00423, 1.99963, 1.99503, 1.99044};
    
    G4MaterialPropertiesTable* reflLayer2_MPT = new G4MaterialPropertiesTable();
    reflLayer2_MPT->AddProperty("REALRINDEX", reflLayer2_Energy, reflLayer2_RIND, false, true);
    reflLayer2_MPT->AddProperty("IMAGINARYRINDEX", reflLayer2_Energy, reflLayer2_KIND, false, true);
    refl2_mat->SetMaterialPropertiesTable(reflLayer2_MPT);
    
    
    // glue
    std::vector<G4double> glue_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
    G4double glue_rindex_const = 1.4303;
    std::vector<G4double> glue_RIND    = { glue_rindex_const, glue_rindex_const, glue_rindex_const, glue_rindex_const, glue_rindex_const};
    G4MaterialPropertiesTable* glue_MPT = new G4MaterialPropertiesTable();
    glue_MPT->AddProperty("RINDEX", glue_Energy, glue_RIND, false, true);
    glue_mat->SetMaterialPropertiesTable(glue_MPT);

    // Air
    std::vector<G4double> air_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
    G4double air_rindex_const = 1.00027717;
    std::vector<G4double> air_RIND    = { air_rindex_const, air_rindex_const, air_rindex_const, air_rindex_const, air_rindex_const};
    G4MaterialPropertiesTable* air_MPT = new G4MaterialPropertiesTable();
    air_MPT->AddProperty("RINDEX", air_Energy, air_RIND, false, true);
    air_mat->SetMaterialPropertiesTable(air_MPT);
    
    // SiPad
    
    std::vector<G4double> SiPad_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
    std::vector<G4double> SiPad_RIND  = { 1.42, 1.42, 1.42, 1.42, 1.42};
    G4MaterialPropertiesTable* SiPad_MPT = new G4MaterialPropertiesTable();
    SiPad_MPT->AddProperty("RINDEX", SiPad_Energy, SiPad_RIND, false, true);
    SiPad_mat->SetMaterialPropertiesTable(SiPad_MPT);
    
    // Borosilicate glass
    
    std::vector<G4double> boro_glass_Energy  = {1.*eV, 1.8246*eV, 1.8657*eV, 1.9088*eV, 1.9539*eV, 2.0011*eV,
        2.0511*eV, 2.1029*eV, 2.1577*eV, 2.2155*eV, 2.2765*eV,
        2.3409*eV, 2.4091*eV, 2.4814*eV, 2.5581*eV, 2.6398*eV,
        2.7268*eV, 2.8198*eV, 2.9193*eV, 3.0261*eV, 3.1410*eV,
        3.265*eV, 10.*eV};
    std::vector<G4double> boro_glass_RIND  = {1.53374486593, 1.53374486593, 1.53152807464, 1.52956876722, 1.52782657556,
        1.52626886013, 1.52486896756, 1.52360493736, 1.52245852827,
        1.52141447577, 1.52045991938, 1.51958395606, 1.51877728882,
        1.51803194746, 1.51734106514, 1.51669869799, 1.51609967864,
        1.51553949635, 1.51501419841, 1.51452030856, 1.51405475908,
        1.51361483398, 1.51361483398};
    G4MaterialPropertiesTable* boro_glass_MPT = new G4MaterialPropertiesTable();
    boro_glass_MPT->AddProperty("RINDEX", boro_glass_Energy, boro_glass_RIND, false, true);
    boro_glass_mat->SetMaterialPropertiesTable(boro_glass_MPT);
    
    // Bialkali Photocathode
    std::vector<G4double> bialkali_Energy  = {1.*eV, 1.8246*eV, 1.8657*eV, 1.9088*eV, 1.9539*eV, 2.0011*eV,
        2.0511*eV, 2.1029*eV, 2.1577*eV, 2.2155*eV, 2.2765*eV,
        2.3409*eV, 2.4091*eV, 2.4814*eV, 2.5581*eV, 2.6398*eV,
        2.7268*eV, 2.8198*eV, 2.9193*eV, 3.0261*eV, 3.1410*eV,
        3.265*eV, 10.*eV};
    std::vector<G4double> bialkali_RIND  = {2.96, 2.96, 2.95, 2.95, 2.95, 2.96, 2.98, 3.01, 3.06, 3.12,
        3.2, 3.26, 3.09, 3, 3, 3, 2.87, 2.7, 2.61, 2.38, 2.18,
        1.92, 1.92};
    
    std::vector<G4double> bialkali_ABSLENGTH = {163.98*nm, 163.98*nm, 155.64*nm, 152.13*nm, 144.38*nm, 133.35*nm, 126.70*nm,
                    111.79*nm, 99.472*nm, 84.082*nm, 68.841*nm, 49.042*nm, 39.031*nm,
                    37.537*nm, 34.770*nm, 27.912*nm, 25.144*nm, 23.343*nm, 22.105*nm,
                    19.080*nm, 18.599*nm, 17.893*nm, 17.893*nm};

    G4MaterialPropertiesTable* bialkali_MPT = new G4MaterialPropertiesTable();
    bialkali_MPT->AddProperty("RINDEX", bialkali_Energy, bialkali_RIND, false, true);
    bialkali_MPT->AddProperty("ABSLENGTH", bialkali_Energy, bialkali_ABSLENGTH, false, true);
    bialkali_mat->SetMaterialPropertiesTable(bialkali_MPT);
                    
    /*
        ##########################################################################################################
        ##########################################################################################################

        TYPE1: CLAIRE
        
        ##########################################################################################################
        ##########################################################################################################

    */

    // -----------> SiPad

    G4int SiPad_copy = 200;

    G4Tubs* solid_SiPad = new G4Tubs("solid_SiPad", 0.*cm, 2.3*cm, 2.5*mm, 0.*deg, 360.*deg);
    G4LogicalVolume* log_SiPad  = new G4LogicalVolume(solid_SiPad, SiPad_mat, "log_SiPad");

    G4VPhysicalVolume* phys_SiPad = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*mm, 22.8*mm), log_SiPad, "phys_SiPad", log_chamber, false, SiPad_copy, true);

    //gm.AddXYZDetector( log_SiPad);
    G4VisAttributes* VisSiPad = new G4VisAttributes(G4Colour::Black());
    log_SiPad->SetVisAttributes(VisSiPad);
    
    // -----------> Glass
    
    G4int Glass_copy = 201;
    
    G4Tubs* solidGlass = new G4Tubs("solidGlass", 0.*cm, 2.3*cm, 2.5*mm, 0.*deg, 360.*deg);
    G4LogicalVolume* logicGlass = new G4LogicalVolume(solidGlass, boro_glass_mat, "logicGlass");
    G4VPhysicalVolume* physGlass = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*mm, 27.8*mm), logicGlass, "Glass", log_chamber, false, Glass_copy, true);
    
    gm.AddXYZDetector( logicGlass);
    G4VisAttributes* VisGlass = new G4VisAttributes(G4Colour::Blue());
    logicGlass->SetVisAttributes(VisGlass);
    
    // -----------> Bialkali
    G4int bialkali_copy = 202;
    
    G4Tubs* solidBialkali = new G4Tubs("solidBialkali", 0.*cm, 2.3*cm, 10*nm, 0.*deg, 360.*deg);
    G4LogicalVolume* logicBialkali = new G4LogicalVolume(solidBialkali, bialkali_mat, "solidBialkali");
    G4VPhysicalVolume* physBialkali = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*mm, 30.30001*mm), logicBialkali, "Bialkali", log_chamber, false, bialkali_copy, true);
    
    gm.AddXYZDetector( logicBialkali);
    G4VisAttributes* VisBialkali = new G4VisAttributes(G4Colour::Green());
    logicBialkali->SetVisAttributes(VisBialkali);

    // -----------> PMT
    G4int PMT_copy = 203;
    
    G4Tubs* solidPMT = new G4Tubs("solidPMT", 0.*cm, 2.3*cm, 1.*mm, 0.*deg, 360.*deg);
    G4LogicalVolume* logicPMT = new G4LogicalVolume(solidPMT, chamber_mat, "logicPMT");
    G4VPhysicalVolume* physPMT = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*mm, 31.30002*mm), logicPMT, "PMT", log_chamber, false, PMT_copy, true);
    
    //gm.AddXYZDetector( logicPMT);
    G4VisAttributes* VisPMT = new G4VisAttributes(G4Colour::Red());
    logicPMT->SetVisAttributes(VisPMT);

    // -----------> AuXDet
    G4int Aux_copy = 204;
    
    G4Tubs* solidAux = new G4Tubs("solidAux", 0.*cm, 2.3*cm, 10.*mm, 0.*deg, 360.*deg);
    G4LogicalVolume* logicAux = new G4LogicalVolume(solidAux, chamber_mat, "logicAux");
    G4VPhysicalVolume* physAux = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*mm, 42.30002*mm), logicAux, "Aux", log_chamber, false, Aux_copy, true);
    
    //gm.AddXYZDetector( logicAux);
    G4VisAttributes* VisAux = new G4VisAttributes(G4Colour::White());
    logicAux->SetVisAttributes(VisAux);
    /*
    // -----------> Air
    G4Tubs* solidAir = new G4Tubs("solidAir", 0.*cm, 2.3*cm, 0.005*mm, 0.*deg, 360.*deg); 
    G4LogicalVolume* logicAir = new G4LogicalVolume(solidAir, air_mat, "logicAir");
    G4VPhysicalVolume* physAir = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*mm, 20.305*mm), logicAir, "Air", log_chamber, false, 300, true);

    // -----------> AirTop
    G4Tubs* solidAirTop = new G4Tubs("solidAirTop", 0.*cm, 2.3*cm, 0.005*mm, 0.*deg, 360.*deg); 
    G4LogicalVolume* logicAirTop = new G4LogicalVolume(solidAirTop, air_mat, "logicAirTop");
    G4VPhysicalVolume* physAirTop = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*mm, 25.315*mm), logicAirTop, "AirTop", log_chamber, false, 301, true);
    */
    
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  -------> CAD volumes

    // CAD path
    G4String cad_path = "";
    gm.config->readInto(cad_path, "GEOM.CAD.PATH");
    G4cout << "GEOM.CAD.PATH: " << cad_path << G4endl;

    
    // -----------> CsI
    auto mesh_CsI = CADMesh::TessellatedMesh::FromPLY(cad_path+"/CsI_block_V3.stl");
    G4VSolid* solid_CsI = mesh_CsI->GetSolid();

    G4int CsI_copy = 100;

    G4LogicalVolume* log_CsI  = new G4LogicalVolume(solid_CsI, csi_mat, "log_CsI");
    
    auto rotation = new G4RotationMatrix();
    rotation->rotateZ(-90*deg);
    rotation->rotateX(-90*deg);
    
    
    G4VPhysicalVolume* phys_CsI = new G4PVPlacement(rotation, G4ThreeVector(46.5*mm, -118.5*mm, 0*mm), log_CsI, "phys_CsI", log_chamber, false, CsI_copy, true);

    //gm.AddXYZDetector( log_CsI);
    G4VisAttributes* VisCsI = new G4VisAttributes(G4Colour::Yellow());
    log_CsI->SetVisAttributes(VisCsI);
    

    // -----------> Teflon
    auto mesh_Teflon = CADMesh::TessellatedMesh::FromPLY(cad_path+"/CsI_block_Teflon.stl");
    G4VSolid* solid_Teflon = mesh_Teflon->GetSolid();

    G4int Teflon_copy = 101;

    G4LogicalVolume* log_Teflon  = new G4LogicalVolume(solid_Teflon, refl1_mat, "log_Teflon");

    G4VPhysicalVolume* phys_Teflon = new G4PVPlacement(rotation, G4ThreeVector(46.5*mm, -118.5*mm, 0*mm), log_Teflon, "phys_Teflon", log_chamber, false, Teflon_copy, true);

    //gm.AddXYZDetector( log_Teflon);
    G4VisAttributes* VisTeflon = new G4VisAttributes(G4Colour::Black());
    log_Teflon->SetVisAttributes(VisTeflon);
    /*
    // -----------> AlMain
    auto mesh_AlMain = CADMesh::TessellatedMesh::FromPLY(cad_path+"/CsI_block_main_structure.stl");
    G4VSolid* solid_AlMain = mesh_AlMain->GetSolid();

    G4int AlMain_copy = 102;

    G4LogicalVolume* log_AlMain  = new G4LogicalVolume(solid_AlMain, refl2_mat, "log_AlMain");

    rotation = new G4RotationMatrix();
    rotation->rotateY(-90*deg);
    rotation->rotateZ(-90*deg);
    G4VPhysicalVolume* phys_AlMain = new G4PVPlacement(rotation, G4ThreeVector(56.5*mm, -125.*mm, -4.*mm), log_AlMain, "phys_AlMain", log_chamber, false, AlMain_copy, true);

    gm.AddXYZDetector( log_AlMain);
    G4VisAttributes* VisAlMain = new G4VisAttributes(G4Colour::Grey());
    log_AlMain->SetVisAttributes(VisAlMain);
    */
    // -----------> AlBot
    auto mesh_AlBot = CADMesh::TessellatedMesh::FromPLY(cad_path+"/CsI_block_structure_bottom.stl");
    G4VSolid* solid_AlBot = mesh_AlBot->GetSolid();

    G4int AlBot_copy = 103;

    G4LogicalVolume* log_AlBot  = new G4LogicalVolume(solid_AlBot, refl2_mat, "log_AlBot");

    G4VPhysicalVolume* phys_AlBot = new G4PVPlacement(rotation, G4ThreeVector(56.5*mm, -125.*mm, -4.*mm), log_AlBot, "phys_AlBot", log_chamber, false, AlBot_copy, true);

    //gm.AddXYZDetector( log_AlBot);
    G4VisAttributes* VisAlBot = new G4VisAttributes(G4Colour::Grey());
    log_AlBot->SetVisAttributes(VisAlBot);
    
    
    /*
    ######### OPTICAL SURFACES ##########
    */
    
    
    // CsI -> Teflon
    
    G4OpticalSurface* panel_OPTSURFACE = new G4OpticalSurface("panel_OPTSURFACE");
    G4LogicalBorderSurface* panel_BORDER = new G4LogicalBorderSurface("panel_BORDER", phys_CsI, phys_Teflon, panel_OPTSURFACE);

    panel_OPTSURFACE->SetModel(unified);
    panel_OPTSURFACE->SetType(dielectric_dielectric);

    G4MaterialPropertiesTable* panel_OPTSURFACE_MPT = new G4MaterialPropertiesTable();
    if (setReflSurfaceType == 1) {
        panel_OPTSURFACE->SetFinish(groundfrontpainted);
        panel_OPTSURFACE->SetSigmaAlpha(2*degree);
        
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    if (setReflSurfaceType == 2) {
        panel_OPTSURFACE->SetFinish(groundfrontpainted);
        panel_OPTSURFACE->SetSigmaAlpha(12*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    if (setReflSurfaceType == 3) {
        panel_OPTSURFACE->SetFinish(groundbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(12*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 1, 1, 1, 1, 1 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND); // refractive index of the air gap
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    if (setReflSurfaceType == 4) {
        panel_OPTSURFACE->SetFinish(polishedfrontpainted);
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        panel_OPTSURFACE->SetSigmaAlpha(1.3*degree); // Janecek at al. 2009
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    if (setReflSurfaceType == 5) {
        panel_OPTSURFACE->SetFinish(groundbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(12*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND);
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    if (setReflSurfaceType == 6) {
        panel_OPTSURFACE->SetFinish(groundbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(40*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND);
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    if (setReflSurfaceType == 7 || setReflSurfaceType == 8) {
        panel_OPTSURFACE->SetFinish(groundbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(40*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 1, 1, 1, 1, 1 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND);
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    /* ***************************************************************************************************************** */
    // fixing groundbackpainted and 100% specular spike, sigma = {2, 12, 40}
    if (setReflSurfaceType == 10){
        panel_OPTSURFACE->SetFinish(groundbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(2*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 1, 1, 1, 1, 1 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND);
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    if (setReflSurfaceType == 11){
        panel_OPTSURFACE->SetFinish(groundbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(12*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 1, 1, 1, 1, 1 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND);
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    if (setReflSurfaceType == 12){
        panel_OPTSURFACE->SetFinish(groundbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(40*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 1, 1, 1, 1, 1 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND);
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    // polishedfrontpainted
    if (setReflSurfaceType == 13) {
        panel_OPTSURFACE->SetFinish(polishedfrontpainted);
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    // polishedbackpainted
    if (setReflSurfaceType == 14) {
        panel_OPTSURFACE->SetFinish(polishedbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(1.3*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 1, 1, 1, 1, 1 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND);
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    // fixing groundbackpainted and sigma = 40 deg, changing refl. to 100% lambertian
    if (setReflSurfaceType == 15){
        panel_OPTSURFACE->SetFinish(groundbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(40*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND);
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    // fixing sigma = 40 deg, changing to groundfrontpainted (refl. lambertian)
    if (setReflSurfaceType == 16){
        panel_OPTSURFACE->SetFinish(groundfrontpainted);
        panel_OPTSURFACE->SetSigmaAlpha(40*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    // polishedbackpainted wit 40 deg
    if (setReflSurfaceType == 17) {
        panel_OPTSURFACE->SetFinish(polishedbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(40*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 1, 1, 1, 1, 1 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND);
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    // fixing polishedbackpainted and sigma = 1.3 deg, changing refl. to 100% lambertian
    if (setReflSurfaceType == 18){
        panel_OPTSURFACE->SetFinish(polishedbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(1.3*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND);
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    // fixing polishedbackpainted and sigma = 40 deg, changing refl. to 100% lambertian
    if (setReflSurfaceType == 19){
        panel_OPTSURFACE->SetFinish(polishedbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(40*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND);
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    if (setReflSurfaceType == 20){
        panel_OPTSURFACE->SetFinish(groundbackpainted);
        panel_OPTSURFACE->SetSigmaAlpha(2*degree); // Janecek at al. 2009
        std::vector<G4double> OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
        std::vector<G4double> OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };
        std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
        std::vector<G4double> OPTSURFACE_RIND    = { 1.00293, 1.00293, 1.00293, 1.00293, 1.00293};
        panel_OPTSURFACE_MPT->AddProperty("RINDEX", OPTSURFACE_Energy, OPTSURFACE_RIND);
        panel_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", OPTSURFACE_Energy, OPTSURFACE_refl);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARLOBECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_Energy, OPTSURFACE_SPECULARSPIKECONSTANT);
        panel_OPTSURFACE_MPT->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_Energy, OPTSURFACE_BACKSCATTERCONSTANT);
        panel_OPTSURFACE_MPT->DumpTable();
    }
    panel_OPTSURFACE->SetMaterialPropertiesTable(panel_OPTSURFACE_MPT);
    
    
    // ------> Teflon-Aluminum foil
    G4OpticalSurface* refl1_OPTSURFACE = new G4OpticalSurface("refl1_OPTSURFACE");
    //G4LogicalBorderSurface* refl1_minusX_BORDER = new G4LogicalBorderSurface("refl1_minusX_BORDER", phys_Teflon, phys_AlMain, refl1_OPTSURFACE);

    refl1_OPTSURFACE->SetType(dielectric_metal);
    refl1_OPTSURFACE->SetFinish(polished);
    refl1_OPTSURFACE->SetModel(unified);
    refl1_OPTSURFACE->SetSigmaAlpha(1.3*degree);
    
    G4MaterialPropertiesTable* refl1_OPTSURFACE_MPT = new G4MaterialPropertiesTable();
    std::vector<G4double> refl1_OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
    //std::vector<G4double> refl1_OPTSURFACE_refl  = { 0.787, 0.787, 0.787, 0.787, 0.787 };
    std::vector<G4double> refl1_OPTSURFACE_refl  = { 0.99, 0.99, 0.99, 0.99, 0.99 };

    refl1_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", refl1_OPTSURFACE_Energy, refl1_OPTSURFACE_refl);
    refl1_OPTSURFACE_MPT->DumpTable();
    
    refl1_OPTSURFACE->SetMaterialPropertiesTable(refl1_OPTSURFACE_MPT);
    
    // ------> to Aux
        
    G4OpticalSurface* OPTSURFACE_EXTERNAL = new G4OpticalSurface("OPTSURFACE_EXTERNAL");
    
    /*
    G4LogicalBorderSurface* SiPad_BORDER_ext = new G4LogicalBorderSurface("SiPad_BORDER_ext", phys_SiPad, phys_chamber, OPTSURFACE_EXTERNAL);
    G4LogicalBorderSurface* Glass_BORDER_ext = new G4LogicalBorderSurface("Glass_BORDER_ext", physGlass, phys_chamber, OPTSURFACE_EXTERNAL);
    G4LogicalBorderSurface* Bialkali_BORDER_ext = new G4LogicalBorderSurface("Bialkali_BORDER_ext", physBialkali, phys_chamber, OPTSURFACE_EXTERNAL);
    G4LogicalBorderSurface* PMT_BORDER_ext = new G4LogicalBorderSurface("PMT_BORDER_ext", physPMT, phys_chamber, OPTSURFACE_EXTERNAL);
    G4LogicalBorderSurface* Aux_BORDER_ext = new G4LogicalBorderSurface("Aux_BORDER_ext", physAux, phys_chamber, OPTSURFACE_EXTERNAL);
    G4LogicalBorderSurface* CsI_BORDER_ext = new G4LogicalBorderSurface("CsI_BORDER_ext", phys_CsI, phys_chamber, OPTSURFACE_EXTERNAL);
    */
    
    G4LogicalBorderSurface* SiPad_Aux_ext = new G4LogicalBorderSurface("SiPad_Aux_ext", phys_SiPad, physAux, OPTSURFACE_EXTERNAL);
    G4LogicalBorderSurface* Glass_Aux_ext = new G4LogicalBorderSurface("Glass_Aux_ext", physGlass, physAux, OPTSURFACE_EXTERNAL);
    G4LogicalBorderSurface* Bialkali_Aux_ext = new G4LogicalBorderSurface("Bialkali_Aux_ext", physBialkali, physAux, OPTSURFACE_EXTERNAL);
    //G4LogicalBorderSurface* PMT_Aux_ext = new G4LogicalBorderSurface("PMT_Aux_ext", physPMT, physAux, OPTSURFACE_EXTERNAL);
    G4LogicalBorderSurface* Teflon_Aux_ext = new G4LogicalBorderSurface("Teflon_Aux_ext", phys_Teflon, physAux, OPTSURFACE_EXTERNAL);
    
    OPTSURFACE_EXTERNAL->SetType(dielectric_metal);
    OPTSURFACE_EXTERNAL->SetFinish(polished);
    OPTSURFACE_EXTERNAL->SetModel(unified);
        
    G4MaterialPropertiesTable* OPTSURFACE_MPT_ext = new G4MaterialPropertiesTable();
    std::vector<G4double> OPTSURFACE_ext_energy = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
    std::vector<G4double> OPTSURFACE_ext_eff  = { 1, 1, 1, 1, 1 };
    std::vector<G4double> OPTSURFACE_ext_refl  = { 0, 0, 0, 0, 0 };

    OPTSURFACE_MPT_ext->AddProperty("EFFICIENCY", OPTSURFACE_ext_energy, OPTSURFACE_ext_eff);
    OPTSURFACE_MPT_ext->AddProperty("REFLECTIVITY", OPTSURFACE_ext_energy, OPTSURFACE_ext_refl);
    OPTSURFACE_MPT_ext->DumpTable();
    OPTSURFACE_EXTERNAL->SetMaterialPropertiesTable(OPTSURFACE_MPT_ext);
    
    // ------> to chamber
    G4OpticalSurface* OPTSURFACE_CHAMBER = new G4OpticalSurface("OPTSURFACE_CHAMBER");
    G4LogicalBorderSurface* SiPad_BORDER_ext = new G4LogicalBorderSurface("SiPad_BORDER_ext", phys_SiPad, phys_chamber, OPTSURFACE_CHAMBER);
    G4LogicalBorderSurface* Glass_BORDER_ext = new G4LogicalBorderSurface("Glass_BORDER_ext", physGlass, phys_chamber, OPTSURFACE_CHAMBER);
    G4LogicalBorderSurface* Bialkali_BORDER_ext = new G4LogicalBorderSurface("Bialkali_BORDER_ext", physBialkali, phys_chamber, OPTSURFACE_CHAMBER);
    G4LogicalBorderSurface* PMT_BORDER_ext = new G4LogicalBorderSurface("PMT_BORDER_ext", physPMT, phys_chamber, OPTSURFACE_CHAMBER);
    G4LogicalBorderSurface* Aux_BORDER_ext = new G4LogicalBorderSurface("Aux_BORDER_ext", physAux, phys_chamber, OPTSURFACE_CHAMBER);
    G4LogicalBorderSurface* CsI_BORDER_ext = new G4LogicalBorderSurface("CsI_BORDER_ext", phys_CsI, phys_chamber, OPTSURFACE_CHAMBER);
    
    OPTSURFACE_CHAMBER->SetType(dielectric_dielectric);
    OPTSURFACE_CHAMBER->SetFinish(groundbackpainted);
    OPTSURFACE_CHAMBER->SetSigmaAlpha(12*degree);
        
    G4MaterialPropertiesTable* OPTSURFACE_MPT_chamber = new G4MaterialPropertiesTable();
    std::vector<G4double> OPTSURFACE_chamber_energy = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
    std::vector<G4double> OPTSURFACE_chamber_refl  = { 1, 1, 1, 1, 1 };
    std::vector<G4double> OPTSURFACE_chamber_rindex  = { 1, 1, 1, 1, 1 };
    std::vector<G4double> OPTSURFACE_SPECULARLOBECONSTANT  = { 0, 0, 0, 0, 0 };
    std::vector<G4double> OPTSURFACE_SPECULARSPIKECONSTANT  = { 1, 1, 1, 1, 1 };
    std::vector<G4double> OPTSURFACE_BACKSCATTERCONSTANT  = { 0, 0, 0, 0, 0 };
    
    OPTSURFACE_MPT_chamber->AddProperty("RINDEX", OPTSURFACE_chamber_energy, OPTSURFACE_chamber_rindex);
    OPTSURFACE_MPT_chamber->AddProperty("REFLECTIVITY", OPTSURFACE_chamber_energy, OPTSURFACE_chamber_refl);
    OPTSURFACE_MPT_chamber->AddProperty("SPECULARLOBECONSTANT", OPTSURFACE_chamber_energy, OPTSURFACE_SPECULARLOBECONSTANT);
    OPTSURFACE_MPT_chamber->AddProperty("SPECULARSPIKECONSTANT", OPTSURFACE_chamber_energy, OPTSURFACE_SPECULARSPIKECONSTANT);
    OPTSURFACE_MPT_chamber->AddProperty("BACKSCATTERCONSTANT", OPTSURFACE_chamber_energy, OPTSURFACE_BACKSCATTERCONSTANT);
    
    OPTSURFACE_MPT_chamber->DumpTable();
    OPTSURFACE_CHAMBER->SetMaterialPropertiesTable(OPTSURFACE_MPT_chamber);
    
    /*
    // ------> to PMT
    
    G4OpticalSurface* pmt_OPTSURFACE = new G4OpticalSurface("pmt_OPTSURFACE");

    G4LogicalBorderSurface* pmt_BORDER = new G4LogicalBorderSurface("pmt_BORDER", physBialkali, physPMT, pmt_OPTSURFACE);

    pmt_OPTSURFACE->SetType(dielectric_metal);
    pmt_OPTSURFACE->SetFinish(polished);
    pmt_OPTSURFACE->SetModel(unified);
    
    G4MaterialPropertiesTable* pmt_OPTSURFACE_MPT = new G4MaterialPropertiesTable();
    std::vector<G4double> pmt_OPTSURFACE_Energy  = { 1.0 * eV, 3.0 * eV, 5.0 * eV, 7.0 * eV, 10 * eV };
    std::vector<G4double> pmt_OPTSURFACE_eff  = { 1, 1, 1, 1, 1 };
    std::vector<G4double> pmt_OPTSURFACE_refl  = { 0, 0, 0, 0, 0 };

    pmt_OPTSURFACE_MPT->AddProperty("EFFICIENCY", pmt_OPTSURFACE_Energy, pmt_OPTSURFACE_eff);
    pmt_OPTSURFACE_MPT->AddProperty("REFLECTIVITY", pmt_OPTSURFACE_Energy, pmt_OPTSURFACE_refl);
    pmt_OPTSURFACE_MPT->DumpTable();
    
    pmt_OPTSURFACE->SetMaterialPropertiesTable(pmt_OPTSURFACE_MPT);
    */
        
    
    return World_phys;
}