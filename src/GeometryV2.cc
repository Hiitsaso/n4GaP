#include <n4-all.hh>

#include "GeometryV2.hh"
#include "n4-utils.hh"

#include "materials.hh"
#include "GeometryV2.hh"
#include "detector.hh"

//~ #include <FTFP_BERT.hh>
//~ #include <G4EmStandardPhysics_option4.hh>
//~ #include <G4LogicalVolume.hh>
//~ #include <G4Material.hh>
//~ #include <G4OpticalPhysics.hh>
//~ #include <G4PVPlacement.hh>
//~ #include <G4RandomDirection.hh>
//~ #include <G4RunManagerFactory.hh>
#include <G4SubtractionSolid.hh>
//~ #include <G4SystemOfUnits.hh>
//~ #include <G4Tubs.hh>
//~ #include <G4UnionSolid.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
//~ #include <G4UserLimits.hh>
//~ #include <G4UserLimits.hh>

//~ #include <G4VSolid.hh>
//~ #include <iostream>
//~ #include <typeinfo>

using vecd = std::vector<G4double>;

const auto world_size = 0.6*m;

G4Material* peek;
G4Material* steel;
G4Material* aluminum;
G4Material* Cu;
G4Material* vacuum;
G4Material* mesh_mat;
G4Material* quartz;
G4Material* tpb;
G4Material* gas;
G4Material* air;
G4Material* teflon;
G4Material* plastic;
G4Material* silicon;
G4Material* Pb;
G4LogicalVolume* world;


field_cage_parameters version2_parameters() {
  field_cage_parameters fcp;
  //PROPERTIES
  fcp.photoe_prob       =   0;
  fcp.pressure          =  10 * bar;
  fcp.temperature       = 293 * kelvin;
  //const auto sc_yield        =  22222./MeV; // Wsc = 45 eV, fr
  fcp.sc_yield          =   1/GeV;
  //~ const auto sc_yield        =  1000./MeV;
  fcp.elifetime         = 1e6 * ms;
  // const auto drift_vel         = 1. * mm/microsecond;
  // const auto drift_transv_diff = 1. * mm/sqrt(cm);
  // const auto drift_long_diff   = .3 * mm/sqrt(cm);
  // const auto el_field          = 16.0 * kilovolt/cm;
  // const auto el_vel            = 3. * mm/microsecond;
  // const auto el_transv_diff    = 1. * mm/sqrt(cm);
  // const auto el_long_diff      = .3 * mm/sqrt(cm);
  
  //MEASUREMENTS
  fcp.vessel_out_rad    = 288./2  *mm;
  fcp.vessel_out_length = 466.79  *mm;
  fcp.vessel_rad        = 276./2  *mm;
  fcp.vessel_length     = 386.471 *mm; 

  fcp.mesh_rad          = 120./2  *mm;
  fcp.mesh_thickn       = 0.075   *mm;
  fcp.mesh_transparency = 0.95;

  fcp.anodeBracket_rad      = 120./2  *mm;
  fcp.anodeBracket_thickn   = 20.     *mm;
  fcp.gateBracket_rad      = 120./2  *mm;
  fcp.gateBracket_thickn   = 30.     *mm;
  fcp.meshBracket_length   = 5.      *mm; //anode=gate
  fcp.cathBracket_rad     = 120./2  *mm;
  fcp.cathBracket_thickn  = 30.     *mm;
  fcp.cathBracket_length  = 4.      *mm; 

  fcp.pmt_rad = 25.4/2  *mm;

  fcp.enclosure_pmt_rad        = 120./2  *mm;
  fcp.enclosure_pmt_thickn     = 8.5     *mm;
  fcp.enclosure_pmt_length     = 113.5   *mm;
  fcp.enclosurevac_pmt_length  = 110.5   *mm;

  fcp.plate_pmt_rad        = 105./2  *mm;
  fcp.plate_pmt_thickn     = 105./2  *mm;
  fcp.plate_pmt_length     = 20      *mm;
  fcp.plateUp_pmt_length   = 15      *mm;
  fcp.plateUp_pmt_thickn   = 21.5    *mm;
  fcp.plateUp_pmt_rad      = 126./2  *mm;
  fcp.plateBottom_pmt_length   = 22      *mm;
  fcp.plateBottom_pmt_rad      = 105./2  *mm;
  
  fcp.pmtHolder_rad        = 115./2  *mm;
  fcp.pmtHolder_length     = 9       *mm;

  fcp.quartz_window_rad    = 108./2  *mm;
  fcp.quartz_window_thickn = 3       *mm;
  fcp.tpb_coating_thickn   = 3       *micrometer;

  fcp.pmt_length = 43.0  *mm;
  
  fcp.TPB_PMTs_rad    = fcp.pmt_rad;
  fcp.TPB_PMTs_length = 3 *micrometer;
  
  fcp.rings_rad    = 165./2  *mm;
  fcp.rings_thickn = 5.      *mm;
  fcp.rings_length = 10.     *mm;
  
  fcp.ring1_rad                   = 40./2   *mm;
  fcp.ring1_thickn                = 19.     *mm;
  fcp.ring2_rad                   = 84./2   *mm;
  fcp.ring2_thickn                = 10.     *mm;
  fcp.ring3_rad                   = 110./2  *mm;
  fcp.ring3_thickn                = 5.      *mm;
  fcp.ring_length                 = 9.      *mm;
  fcp.ring_to_ring                = 3.      *mm;
  fcp.ring_bottom_to_ring0_bottom = 5.      *mm;
  fcp.ring0_length                = 5.      *mm;
  
  fcp.teflon_cage_rad    = 120./2  *mm;
  fcp.teflon_cage_thickn = 5.      *mm;
  fcp.teflon_cage_length = 87.     *mm;
  
  fcp.TPB_tefloncage_length = fcp.teflon_cage_length - (fcp.ring_length - fcp.cathBracket_length);
  fcp.TPB_tefloncage_thickn = 3 *mm;
  fcp.TPB_tefloncage_rad    = fcp.teflon_cage_rad;
  
  fcp.encapsulation_rad    = 38./2  *mm;
  fcp.encapsulation_length = 2   *mm; 
  
  fcp.SiPMs_long         = 6.   *mm; 
  fcp.SiPMs_short        = 6.   *mm; 
  fcp.SiPMs_thickn       = 2.   *mm; //¿?
  fcp.SiPMs_cage_long    = 15.  *mm; 
  fcp.SiPMs_cage_short   = 15.  *mm; 
  fcp.SiPMs_cage_thickn  = 2.   *mm; 
  fcp.anode_to_SiPMs     = 5.   *mm; 
  fcp.SiPMs_ef_xy        = 6.   *mm;
  fcp.SiPMs_ef_between   = 0.5  *mm;
  
  fcp.SiPM_between_long  = 5. *mm;
  fcp.SiPM_between_short = 5. *mm;
  fcp.SiPM_number        = 4;
  
  fcp.Pb_box_lengt_xy = 100 *mm;
  fcp.Pb_box_lengt_z  = 50  *mm;
  fcp.Pb_box_rel_pos  = 2   *cm;
  
  //S1 AND S2 LENGTHS
  fcp.drift_length = 87 *mm; // ||   |-|-------------------||
  fcp.el_length    = 10 *mm; // ||---| |                   ||
  
  fcp.S2_rad    = fcp.anodeBracket_rad + fcp.anodeBracket_thickn;
  fcp.S2_lenght = fcp.el_length - 2*fcp.mesh_thickn; //Since the meshes are outside of the anode/gate Brackets
  
  fcp.S1_rad    = fcp.teflon_cage_rad - fcp.TPB_tefloncage_thickn; 
  fcp.S1_lenght = fcp.drift_length;  
  
  
  //POSITIONS(faltan por medir)
  fcp.vessel_z = fcp.vessel_length/2 - (163.5*mm + fcp.meshBracket_length);
  
  fcp.anodeBracket_z  = -fcp.vessel_z; //From the vessel (and at the origin or the world)
  fcp.anode_z      = fcp.anodeBracket_z + fcp.meshBracket_length/2 + fcp.mesh_thickn/2; 
  fcp.gateBracket_z  = fcp.anodeBracket_z + fcp.meshBracket_length   + fcp.el_length;
  fcp.gate_z         = fcp.gateBracket_z - fcp.meshBracket_length/2 - fcp.mesh_thickn/2;
  fcp.cathBracket_z = fcp.gateBracket_z + fcp.meshBracket_length/2 + fcp.drift_length + fcp.cathBracket_length/2;
  
  fcp.teflon_cage_z    = fcp.gateBracket_z  +  fcp.meshBracket_length/2 + fcp.teflon_cage_length/2; 
  fcp.TPB_tefloncage_z = fcp.teflon_cage_z  - (fcp.ring_length - fcp.cathBracket_length)/2;
  fcp.long_ring_z      = fcp.cathBracket_z - (fcp.ring_length - fcp.cathBracket_length)/2; 
  fcp.teflon_ring_z    = fcp.cathBracket_z - (fcp.ring_length - fcp.cathBracket_length)/2 - fcp.ring_length/4; 
  fcp.cath_ring_z   = fcp.cathBracket_z - (fcp.ring_length - fcp.cathBracket_length)/2 + fcp.ring_length/4; 
  fcp.ring0_z          = fcp.long_ring_z    + fcp.rings_length/2  -  fcp.ring_bottom_to_ring0_bottom - fcp.ring0_length/2; 
  
  //~ fcp.pmt_z             = - (12.*mm + fcp.meshBracket_length/2 + fcp.pmt_length/2) + fcp.anodeBracket_z; //old position of the PMTs
  fcp.pmt_z                 = - (2.*mm + fcp.meshBracket_length/2 + fcp.pmt_length/2) + fcp.anodeBracket_z;  //new position of the PMTs (2mm from the cath)
  fcp.plateUp_pmt_z         = - ( 150.5*mm + fcp.meshBracket_length/2 + fcp.plateUp_pmt_length/2) + fcp.anodeBracket_z;
  fcp.enclosure_pmt_z       = fcp.plateUp_pmt_z + fcp.enclosure_pmt_length/2 + fcp.plateUp_pmt_length/2;
  fcp.plate_pmt_z           = fcp.enclosure_pmt_z + fcp.enclosure_pmt_length/2 + fcp.plate_pmt_length/2;
  fcp.PMTplateBottom1_pos_z = fcp.pmt_z;
  
  fcp.SiPMs_z = fcp.anodeBracket_z - fcp.meshBracket_length/2 - fcp.anode_to_SiPMs;
   
  fcp.drift_z = fcp.gateBracket_z + fcp.meshBracket_length + fcp.drift_length/2; //not useful
  fcp.el_z    = fcp.anodeBracket_z + fcp.meshBracket_length + fcp.el_length/2;    //not useful
  
  fcp.S2_z    = fcp.anode_z + fcp.mesh_thickn/2 + fcp.S2_lenght/2; 
  fcp.S1_z    = fcp.gate_z + fcp.mesh_thickn/2 + fcp.S1_lenght/2;
  
  fcp.encapsulation_z = fcp.ring0_z + (fcp.ring0_length + fcp.encapsulation_length)/2;
  
  fcp.Pb_box_z = fcp.vessel_z + fcp.vessel_length/2 - fcp.Pb_box_rel_pos - fcp.Pb_box_lengt_xy/2;

  return fcp;
}

void ensure_initialized(field_cage_parameters const & fcp) {
  static bool initialized = false;
  if (initialized) { return; }
  initialized = true;

  Cu     = n4::material("G4_Cu"); 
  vacuum = n4::material("G4_Galactic");
  steel  = steel_with_properties();
  aluminum = aluminum_with_properties();
  gas      = GAr_with_properties(fcp.pressure, fcp.temperature, fcp.sc_yield, fcp.elifetime);
  mesh_mat = FakeDielectric_with_properties(gas, "mesh_mat",
                                            fcp.pressure, fcp.temperature, fcp.mesh_transparency, fcp.mesh_thickn,
                                            fcp.sc_yield, fcp.elifetime, fcp.photoe_prob);
  peek    = peek_with_properties();
  quartz  = quartz_with_properties();
  tpb     = TPB_with_properties();
  teflon  = teflon_with_properties();
  plastic = plastic_with_properties();
  silicon = silicon_with_properties();
  Pb      = Pb_with_properties();
  air     = air_with_properties();
  world = n4::box("world").cube(world_size).volume(vacuum);
}

G4LogicalVolume* get_world(field_cage_parameters const & fcp) {
  ensure_initialized(fcp);
  return world;
}

void TEFLON_place_optical_surface_between(G4PVPlacement* one, G4PVPlacement* two, G4String name){
    static G4OpticalSurface* optical_surface = nullptr;
    if (! optical_surface) {
        optical_surface = new G4OpticalSurface(name);
        // Values from same paper as above ("Optimization of Parameters...")
        // "groundfrontpainted" (I think) only considers whether the photon is reflected or absorbed, so there will be no shine through visible in the simulation
        // According to the docs, for UNIFIED, dielectric_dielectric surfaces only the Lambertian reflection is turned on
        optical_surface -> SetType(dielectric_dielectric);
        optical_surface -> SetModel(unified);
        optical_surface -> SetFinish(polished);
        optical_surface -> SetSigmaAlpha(0.0);

        vecd pp = n4::scale_by( eV, {8.8560, 8.2656, 7.7490, 7.2932, 6.8880, 5.6356, 5.1660, 4.7686, 4.4280, 4.1328, 3.5424, 3.0996, 2.4797});
        vecd R  =                   {0.05, 0.05, 0.15, 0.30, 0.60, .88, .90, .93, .95, .98, .98, .98, .98};
        optical_surface -> SetMaterialPropertiesTable(
            n4::material_properties{}
            .add("REFLECTIVITY", pp, R)
            .done());
    }
    new G4LogicalBorderSurface(name, one, two, optical_surface);
}

void place_skin_optical_surface_between(G4LogicalVolume* one, G4String name){
    //Para que no se queden atrapados en el tpb
    static G4OpticalSurface* optical_surface = nullptr;
    if (! optical_surface) {
        optical_surface = new G4OpticalSurface(name);
        optical_surface -> SetType(dielectric_dielectric);
        optical_surface -> SetModel(glisur);
        optical_surface -> SetFinish(ground);
        optical_surface -> SetSigmaAlpha(0.01);
    }
    
    vecd pp = n4::scale_by( eV, {2.4797, 2.8178, 2.9520, 3.0996, 3.2627, 3.3509, 3.4440, 3.6466, 3.8745});
    vecd R  =                   {.84, .82, .82, .82, .52, .36, .26, .30, .32};
    vecd T  =                   {.16, .17, .18, .19, .16, .14, .13, .11, .12};

    optical_surface -> SetMaterialPropertiesTable(
            n4::material_properties{}
            .add("REFLECTIVITY",  pp, R)
            .add("TRANSMITTANCE", pp, T)
            .done());
            
    new G4LogicalSkinSurface(name, one, optical_surface);
}

void TPB_place_optical_surface_between(G4PVPlacement* one, G4PVPlacement* two, G4String name){
    static G4OpticalSurface* optical_surface = nullptr;
    if (! optical_surface) {
        optical_surface = new G4OpticalSurface(name);
        optical_surface -> SetType(dielectric_dielectric);
        optical_surface -> SetModel(glisur);
        optical_surface -> SetFinish(ground);
        optical_surface -> SetSigmaAlpha(0.01);
    }
    
    vecd pp = n4::scale_by( eV, {2.4797, 2.8178, 2.9520, 3.0996, 3.2627, 3.3509, 3.4440, 3.6466, 3.8745});
    vecd R  =                   {.84, .82, .82, .82, .52, .36, .26, .30, .32};
    vecd T  =                   {.16, .17, .18, .19, .16, .14, .13, .11, .12};

    optical_surface -> SetMaterialPropertiesTable(
            n4::material_properties{}
            .add("REFLECTIVITY",  pp, R)
            .add("TRANSMITTANCE", pp, T)
            .done());
            
    new G4LogicalBorderSurface(name, one, two, optical_surface);
}

void place_pmt_holder_in(G4LogicalVolume* vessel, G4PVPlacement* vessel_placement, field_cage_parameters const & fcp, bool TPBon, G4String detector) {     
  //Upper steel plate at the pmt clad
  n4::tubs("PMTplateUp").r_inner(fcp.plateUp_pmt_rad).r_delta(fcp.plateUp_pmt_thickn).z(fcp.plateUp_pmt_length).place(steel).in(vessel).at_z(fcp.plateUp_pmt_z).check_overlaps().now();
  
  //PMT clad
  n4::tubs("EnclosurePMT").r_inner(fcp.enclosure_pmt_rad).r_delta(fcp.enclosure_pmt_thickn).z(fcp.enclosure_pmt_length).place(steel).in(vessel).at_z(fcp.enclosure_pmt_z).check_overlaps().now();
  
  // Steel plate attached where the peek holders are attached
  n4::tubs("PMTplateBottom0").r_inner(fcp.plate_pmt_rad).r_delta(fcp.plate_pmt_thickn).z(fcp.plate_pmt_length).place(steel).in(vessel).at_z(fcp.plate_pmt_z).check_overlaps().now();
  
  //Build PMT and PMTplateBottom1
  auto solid_pmt = n4::tubs("solid_PMT").r(fcp.pmt_rad).z(fcp.pmt_length).solid(); //Hamamatsu pmt length: 43*mm | STEP pmt gap length: 57.5*mm
  G4LogicalVolume* logic_pmt = new G4LogicalVolume(solid_pmt, aluminum, "PMT");
  
  // Position pairs (x,Y) for PMTs
  std::vector <float> pmt_PsX={-15.573,  20.68 , -36.253, 0, 36.253, -20.68 , 15.573};
  std::vector <float> pmt_PsY={-32.871, -29.922,  -2.949, 0,  2.949,  29.922, 32.871};
 
  G4ThreeVector pos_pmt = {0, 0, 0};
  G4ThreeVector pos_tpb = {0, 0, 0};
  
  std::vector<std::vector<G4ThreeVector>> all_pos_pmt;
  
  for (G4int i = 0; i < G4int(pmt_PsX.size()); i++) {
    pos_pmt = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm, fcp.pmt_z);
    
    G4ThreeVector temp_pos_pmt;
    std::vector<G4ThreeVector> temp_vector; 
    temp_pos_pmt = {pos_pmt[0], pos_pmt[1], 0.};
    temp_vector.push_back(temp_pos_pmt); 
    all_pos_pmt.push_back(temp_vector);   
    
    if (detector == "PMT"){
      G4PVPlacement* PMTs = n4::place(logic_pmt).in(vessel).at(pos_pmt).copy_no(i).check_overlaps().now();
      //TPB on top of the PMTs
      if (TPBon == true){
	    pos_tpb = {pos_pmt[0], pos_pmt[1], pos_pmt[2] + fcp.pmt_length/2 + fcp.TPB_PMTs_length/2};
	    auto TPB_in_PMTs = n4::tubs("TPBinPMTs").r(fcp.TPB_PMTs_rad).z(fcp.TPB_PMTs_length).volume(tpb);
	    auto TPB_in_PMTs_placement = n4::place(TPB_in_PMTs).in(vessel).at(pos_tpb).copy_no(i).check_overlaps().now();
	    
        TPB_place_optical_surface_between(TPB_in_PMTs_placement, vessel_placement, "OpticalSurfaceTPBtoVessel");
        TPB_place_optical_surface_between(TPB_in_PMTs_placement, PMTs, "OpticalSurfaceTPBtoPMTs");
	  }
    }		
  }

  auto PMtplateBottom1 =   
  /*      */n4::tubs("solid_PMTplateBottom1").r(fcp.plateBottom_pmt_rad).z(fcp.plateBottom_pmt_length)
  .subtract(solid_pmt).at(all_pos_pmt[0][0])
  .subtract(solid_pmt).at(all_pos_pmt[1][0])
  .subtract(solid_pmt).at(all_pos_pmt[2][0])
  .subtract(solid_pmt).at(all_pos_pmt[3][0])
  .subtract(solid_pmt).at(all_pos_pmt[4][0])
  .subtract(solid_pmt).at(all_pos_pmt[5][0])
  .subtract(solid_pmt).at(all_pos_pmt[6][0]) 
  .name("PMTplateBottom1").place(steel).in(vessel).at_z(fcp.PMTplateBottom1_pos_z).check_overlaps().now();
  
  if (detector == "SiPM"){   
    G4int N = fcp.SiPM_number;
    auto max = 0.;
  
    //~ if (N%2 == 0){ 
	 
	  //~ N = fcp.SiPM_number - 1;  //par
	  //~ G4int lim = (2*fcp.SiPM_number - 5)/4;
  
      //~ for (G4int i = -lim/2; i < lim/2 + 1; ++i) {
		//~ G4double posX = (fcp.SiPMs_cage_long + fcp.SiPM_between_long)*fcp.SiPM_number/2;  //fijo
		//~ G4double posY = fcp.SiPMs_cage_short*i + fcp.SiPM_between_short*i;
		//~ G4ThreeVector pos  = {posX, posY, fcp.SiPMs_z};
		//~ G4ThreeVector pos_ = {-posX, posY, fcp.SiPMs_z};
	    //~ n4::place(logic_4SiPM).in(vessel).at(pos).check_overlaps().now();
	    //~ n4::place(logic_4SiPM).in(vessel).at(pos_).check_overlaps().now();
	  //~ }
      //~ for (G4int j = -lim/2; j < lim/2 + 1; ++j) {
		//~ G4double posY = (fcp.SiPMs_cage_short + fcp.SiPM_between_short)*fcp.SiPM_number/2;  //fijo
		//~ G4double posX =  fcp.SiPMs_cage_long*j + fcp.SiPM_between_long*j;
		//~ G4ThreeVector pos = {posX, posY, fcp.SiPMs_z};
		//~ G4ThreeVector pos_ = {posX, -posY, fcp.SiPMs_z};
	    //~ n4::place(logic_4SiPM).in(vessel).at(pos).check_overlaps().now();
	    //~ n4::place(logic_4SiPM).in(vessel).at(pos_).check_overlaps().now();
	  //~ }
    //~ }  
    
    G4LogicalVolume* logic_4SiPM = n4::box("4SiPM").x(fcp.SiPMs_cage_long).y(fcp.SiPMs_cage_short).z(fcp.SiPMs_thickn).volume(plastic);
	G4Box* solid_4SiPM_ef = n4::box("4SiPMef").x(fcp.SiPMs_ef_xy).y(fcp.SiPMs_ef_xy).z(fcp.SiPMs_thickn).solid();
    G4LogicalVolume* logic_4SiPM_ef = n4::box("4SiPMef").x(fcp.SiPMs_ef_xy).y(fcp.SiPMs_ef_xy).z(fcp.SiPMs_thickn)
									  .add(solid_4SiPM_ef).at({-(fcp.SiPMs_ef_xy + fcp.SiPMs_ef_between), -(fcp.SiPMs_ef_xy + fcp.SiPMs_ef_between), 0.}) 
									  .add(solid_4SiPM_ef).at({-(fcp.SiPMs_ef_xy + fcp.SiPMs_ef_between),                                        0., 0.}) 
									  .add(solid_4SiPM_ef).at({                                       0., -(fcp.SiPMs_ef_xy + fcp.SiPMs_ef_between), 0.}) 
									  .volume(silicon);
	auto sensitive_detector = new n4::sensitive_detector("4Detector", process_hits);
    logic_4SiPM_ef -> SetSensitiveDetector(sensitive_detector); 
    
    for (G4int i = -N/2; i < N/2 + 1; ++i) {
      for (G4int j = -N/2; j < N/2 + 1; ++j) {
		if (i != 0 and j!=0){
			auto i_new = 0.;
			auto j_new = 0.;
			if (i>=0 and j>=0){
				i_new = i - 0.5;
				j_new = j - 0.5;
			}else if (i<=0 and j>=0){
				i_new = i + 0.5;
				j_new = j - 0.5;				
			}else if (i<=0 and j<=0){
				i_new = i + 0.5;
				j_new = j + 0.5;
			}else{
				i_new = i - 0.5;
				j_new = j + 0.5;
			}
			G4double posX = fcp.SiPMs_cage_long*i_new + fcp.SiPM_between_long*i_new;
			G4double posY = fcp.SiPMs_cage_short*j_new + fcp.SiPM_between_short*j_new;
			G4ThreeVector pos = {posX, posY, fcp.SiPMs_z};
			n4::place(logic_4SiPM).in(vessel).at(pos).check_overlaps().now();
			n4::place(logic_4SiPM_ef).in(logic_4SiPM).at({(fcp.SiPMs_ef_xy + fcp.SiPMs_ef_between)/2, (fcp.SiPMs_ef_xy + fcp.SiPMs_ef_between)/2, 0.}).check_overlaps().now();

		}
      }
   } 
   
    //~ G4LogicalVolume* logic_OneSiPM = n4::box("OneSiPM").xy(2*fcp.mesh_rad).z(fcp.SiPMs_thickn).volume(silicon);
    //~ auto sensitive_detector = new n4::sensitive_detector("OneDetector", process_hits);
    //~ logic_OneSiPM -> SetSensitiveDetector(sensitive_detector);
    //~ n4::place(logic_OneSiPM).in(vessel).at_z(fcp.SiPMs_z).check_overlaps().now();
    
  }
}

void place_cage_in(G4LogicalVolume* vessel, field_cage_parameters const & fcp) { 
  //Anode bracket
  n4::tubs("AnodeBracket").r_inner(fcp.anodeBracket_rad).r_delta(fcp.anodeBracket_thickn).z(fcp.meshBracket_length).place(steel).in(vessel).at_z(fcp.anodeBracket_z).check_overlaps().now();
  
  //Anode mesh 
  n4::tubs("Anode").r(fcp.mesh_rad).z(fcp.mesh_thickn).place(mesh_mat).in(vessel).at_z(fcp.anode_z).check_overlaps().now();

  //"Cath" rings (only half of the rings plus the one inside)
  n4::tubs("CathRing0").r(      fcp.ring1_rad).z(fcp.ring0_length).place(steel).in(vessel).at_z(fcp.ring0_z).check_overlaps().now();
  n4::tubs("CathRing1").r_inner(fcp.ring1_rad).r_delta(fcp.ring1_thickn).z(fcp.ring_length/2).place(steel).in(vessel).at_z(fcp.cath_ring_z).check_overlaps().now();
  n4::tubs("CathRing2").r_inner(fcp.ring2_rad).r_delta(fcp.ring2_thickn).z(fcp.ring_length/2).place(steel).in(vessel).at_z(fcp.cath_ring_z).check_overlaps().now();
  n4::tubs("CathRing3").r_inner(fcp.ring3_rad).r_delta(fcp.ring3_thickn).z(fcp.ring_length/2).place(steel).in(vessel).at_z(fcp.cath_ring_z).check_overlaps().now();
  
  //Gate bracket
  n4::tubs("GateBracket").r_inner(fcp.gateBracket_rad).r_delta(fcp.gateBracket_thickn).z(fcp.meshBracket_length).place(steel).in(vessel).at_z(fcp.gateBracket_z).check_overlaps().now();
  
  //Gate mesh
  n4::tubs("Gate").r(fcp.mesh_rad).z(fcp.mesh_thickn).place(mesh_mat).in(vessel).at_z(fcp.gate_z).check_overlaps().now();

  //Cath bracket
  n4::tubs("CathBracket").r_inner(fcp.cathBracket_rad).r_delta(fcp.cathBracket_thickn).z(fcp.cathBracket_length).place(steel).in(vessel).at_z(fcp.cathBracket_z).check_overlaps().now();
}

void place_teflon_cage_in(G4LogicalVolume* vessel, G4PVPlacement* vessel_placement, field_cage_parameters const & fcp, bool TPBon, bool OpticalSurfaceON) { 
  //Teflon cage
  G4PVPlacement* TeflonCage = n4::tubs("TeflonCage").r_inner(fcp.teflon_cage_rad).r_delta(fcp.teflon_cage_thickn).z(fcp.teflon_cage_length).place(teflon).in(vessel).at_z(fcp.teflon_cage_z).check_overlaps().now();  
  //Teflon rings (only half of the rings)
  G4PVPlacement* Teflon_ring_1 = n4::tubs("TeflonRing1").r_inner(fcp.ring1_rad).r_delta(fcp.ring1_thickn).z(fcp.ring_length/2).place(teflon).in(vessel).at_z(fcp.teflon_ring_z).check_overlaps().now();
  G4PVPlacement* Teflon_ring_2 = n4::tubs("TeflonRing2").r_inner(fcp.ring2_rad).r_delta(fcp.ring2_thickn).z(fcp.ring_length/2).place(teflon).in(vessel).at_z(fcp.teflon_ring_z).check_overlaps().now();
  G4PVPlacement* Teflon_ring_3 = n4::tubs("TeflonRing3").r_inner(fcp.ring3_rad).r_delta(fcp.ring3_thickn).z(fcp.ring_length/2).place(teflon).in(vessel).at_z(fcp.teflon_ring_z).check_overlaps().now();
  
  //TPB inside of the teflon cage
  if (TPBon == true){
    G4PVPlacement* TPB_in_TeflonCage = n4::tubs("TPBinTefloncage").r(fcp.TPB_tefloncage_rad).r_delta(fcp.TPB_tefloncage_thickn).z(fcp.TPB_tefloncage_length).place(tpb).in(vessel).at_z(fcp.TPB_tefloncage_z).check_overlaps().now();    
	  if (OpticalSurfaceON == true){
        TEFLON_place_optical_surface_between(TPB_in_TeflonCage,    TeflonCage, "OpticalSurfaceTPBToTeflonCage");
        TEFLON_place_optical_surface_between( vessel_placement, Teflon_ring_1, "OpticalSurfaceVesselToTeflonRing1");
        TEFLON_place_optical_surface_between( vessel_placement, Teflon_ring_2, "OpticalSurfaceVesselToTeflonRing2");
        TEFLON_place_optical_surface_between( vessel_placement, Teflon_ring_3, "OpticalSurfaceVesselToTeflonRing3");
        
        TPB_place_optical_surface_between(TPB_in_TeflonCage,    vessel_placement, "OpticalSurfaceTPBToVessel");
      }
  //~ } else if (TPBon == false && OpticalSurfaceON == true){
	  //~ TEFLON_place_optical_surface_between(vessel_placement,    TeflonCage, "OpticalSurfaceVesseltoTeflonCage");	  
	  //~ TEFLON_place_optical_surface_between(vessel_placement, Teflon_ring_1,   "OpticalSurfaceTPBtoTeflonRing1");	  
	  //~ TEFLON_place_optical_surface_between(vessel_placement, Teflon_ring_2,   "OpticalSurfaceTPBtoTeflonRing2");	  
	  //~ TEFLON_place_optical_surface_between(vessel_placement, Teflon_ring_3,   "OpticalSurfaceTPBtoTeflonRing3");	  
  }
}

void place_encapsulation_in(G4LogicalVolume* vessel, field_cage_parameters const & fcp) {
 n4::tubs("RadioactiveSourceEncapsulation").r(fcp.encapsulation_rad).z(fcp.encapsulation_length).place(plastic).in(vessel).at_z(fcp.encapsulation_z).check_overlaps().now();
 //~ n4::tubs("RadioactiveSourceEncapsulation").r(fcp.encapsulation_rad).r_delta(fcp.encapsulation_thickn).z(fcp.encapsulation_length).place(plastic).in(vessel).at_z(fcp.encapsulation_z).check_overlaps().now();
}

void place_rings_in(G4LogicalVolume* vessel, field_cage_parameters const & fcp) { 
  auto logic_ring = n4::tubs("ring").r_inner(fcp.rings_rad).r_delta(fcp.rings_thickn).z(fcp.rings_length).volume(Cu);
  
  auto gate_to_ring = 5. *mm;
  auto ring_to_ring = 3. *mm;
  auto ring_z       = 0. *mm;
  
  for (G4int i = 0; i < 6; i++) {
    ring_z = fcp.gateBracket_z + fcp.meshBracket_length + gate_to_ring + fcp.rings_length/2*(2*i + 1) + ring_to_ring*i; 
    n4::place(logic_ring).in(vessel).at_z(ring_z).copy_no(i).check_overlaps().now();
  }
}

void place_S1_and_S2_in(G4LogicalVolume* vessel, field_cage_parameters const & fcp){
  //S2 zone
  n4::tubs("S2").r(fcp.S2_rad).z(fcp.S2_lenght).place(gas).in(vessel).at_z(fcp.S2_z).check_overlaps().now();
	
  //S1 zone
  auto S1_logic = 	  
  /*      */n4::tubs("S1_full").r(fcp.S1_rad).z(fcp.S1_lenght)
  .subtract(n4::tubs("CathRing_0_tosubtact").r(fcp.ring1_rad).z(fcp.ring0_length)).at_z(fcp.ring0_z - fcp.S1_z) //Relative positions
  .name("S1").volume(gas);
  
  //~ auto sensitive_detector_2 = new n4::sensitive_detector("S1Detector", process_hits_genratorCHECK);
  //~ S1_logic -> SetSensitiveDetector(sensitive_detector_2);
   
  n4::place(S1_logic).in(vessel).at_z(fcp.S1_z).check_overlaps().now();
}

void place_Pb_box_in(G4LogicalVolume* vessel, field_cage_parameters const & fcp){
  n4::box("PbBox").xy(fcp.Pb_box_lengt_xy).z(fcp.Pb_box_lengt_z).place(Pb).in(vessel).at_z(fcp.Pb_box_z).check_overlaps().now();
}

G4PVPlacement* GeometryV2() {
  field_cage_parameters fcp = version2_parameters();
  ensure_initialized(fcp);
  
  auto vessel_steel =   
  /*      */n4::tubs("vessel_steel_outer").r(fcp.vessel_out_rad).z(fcp.vessel_out_length)
  .subtract(n4::tubs("vessel_steel_inner").r(    fcp.vessel_rad).z(fcp.vessel_length))
  .name("vessel_steel").place(steel).in(world).at_z(fcp.vessel_z).check_overlaps().now();
      
  auto vessel           = n4::tubs("GasVessel").r(fcp.vessel_rad).z(fcp.vessel_length).volume(gas);
  auto vessel_placement = n4::place(vessel).in(world).at_z(fcp.vessel_z).check_overlaps().now();
  
  //bool1 = TPBon ; bool2 = OpticalSurfaceON
  place_pmt_holder_in(vessel, vessel_placement, fcp, true, "SiPM"); //SiPM or PMT
  //~ place_cage_in(vessel, fcp);
  //~ place_teflon_cage_in(vessel, vessel_placement, fcp, true, true);
  //~ place_rings_in(vessel, fcp);
  //place_encapsulation_in(vessel, fcp);
  //~ place_S1_and_S2_in(vessel, fcp);
  //place_Pb_box_in(vessel, fcp);
  
  return n4::place(world).now();
}

G4PVPlacement* GeometryV2_TEST() {
  field_cage_parameters fcp = version2_parameters();
  ensure_initialized(fcp);
  
  //TPB análisis
  //~ auto worldbox = n4::box("Box").xy(world_size).z(world_size).place(gas).in(world).at_z(0.).check_overlaps().now();
  //~ auto box2     = n4::box("Box2").xy(world_size).z(3. *micrometer).volume(tpb);
  //~ // auto box2     = n4::box("Box2").xy(world_size).z(3. *mm).volume(tpb);
  //~ auto box2_plc = n4::place(box2).in(worldbox).at_z(0.).check_overlaps().now();
  //~ // place_skin_optical_surface_between(box2, "opticalskintest");
  //~ // TPB_place_optical_surface_between(worldbox, box2_plc, "opticalsurfacetest1");
  //~ TPB_place_optical_surface_between(box2_plc, worldbox, "opticalsurfacetest2");
  
  //Teflón análisis
  //TPB análisis
  auto worldbox = n4::box("Box").xy(world_size).z(world_size).place(gas).in(world).at_z(0.).check_overlaps().now();
  auto box2     = n4::box("Box2").xy(world_size).z(3. *micrometer).volume(teflon);
  auto box3     = n4::box("Box2").xy(world_size).z(world_size/4).volume(steel);
  //~ auto box2     = n4::box("Box2").xy(world_size).z(3. *mm).volume(tpb);
  auto box2_plc = n4::place(box2).in(worldbox).at_z(0.).check_overlaps().now();
  auto box3_plc = n4::place(box3).in(worldbox).at_z(-world_size/8).check_overlaps().now();
  //~ place_skin_optical_surface_between(box2, "opticalskintest");
  TEFLON_place_optical_surface_between(worldbox, box2_plc, "opticalsurfacetest1");
  //~ place_optical_surface_between(box2_plc, worldbox, "opticalsurfacetest2");
  
  return n4::place(world).now();
}

