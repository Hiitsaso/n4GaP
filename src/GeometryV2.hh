#ifndef GEOMETRYV2_H
#define GEOMETRYV2_H

struct field_cage_parameters {
  //PROPERTIES
  G4double photoe_prob;
  G4double pressure;
  G4double temperature;
  G4double sc_yield;
  G4double elifetime;
  
  //MEASUREMENTS
  G4double vessel_out_rad;
  G4double vessel_out_length;
  G4double vessel_rad;
  G4double vessel_length;
  
  G4double mesh_rad;
  G4double mesh_thickn;
  G4double mesh_transparency;
    
  G4double anodeBracket_rad;
  G4double anodeBracket_thickn;
  G4double gateBracket_rad;
  G4double gateBracket_thickn;
  G4double meshBracket_length;
  G4double cathBracket_rad;
  G4double cathBracket_thickn;
  G4double cathBracket_length;
  
  G4double pmt_rad;
  
  G4double enclosure_pmt_rad;
  G4double enclosure_pmt_thickn;
  G4double enclosure_pmt_length;
  G4double enclosurevac_pmt_length;
  
  G4double plate_pmt_rad;
  G4double plate_pmt_thickn;
  G4double plate_pmt_length;
  G4double plateUp_pmt_length;
  G4double plateUp_pmt_thickn;
  G4double plateUp_pmt_rad;
  G4double plateBottom_pmt_length;
  G4double plateBottom_pmt_rad;
  
  G4double pmtHolder_rad;
  G4double pmtHolder_length;
  
  G4double quartz_window_rad;
  G4double quartz_window_thickn;
  G4double tpb_coating_thickn;
  
  G4double pmt_length;
  
  G4double TPB_PMTs_rad;
  G4double TPB_PMTs_length;
  
  G4double rings_rad;
  G4double rings_thickn;
  G4double rings_length;
  
  G4double ring1_rad;
  G4double ring1_thickn;
  G4double ring2_rad;
  G4double ring2_thickn;
  G4double ring3_rad;
  G4double ring3_thickn;
  G4double ring_length;
  G4double ring_to_ring;
  G4double ring_bottom_to_ring0_bottom;
  G4double ring0_length;
  
  G4double teflon_cage_rad;
  G4double teflon_cage_thickn;
  G4double teflon_cage_length;
  
  G4double TPB_tefloncage_rad;
  G4double TPB_tefloncage_length;
  G4double TPB_tefloncage_thickn;
  
  G4double encapsulation_rad;
  G4double encapsulation_length;
  G4double encapsulation_thickn;
  
  G4double SiPMs_long;
  G4double SiPMs_short;
  G4double SiPMs_thickn;
  G4double SiPMs_cage_long;
  G4double SiPMs_cage_short;
  G4double SiPMs_cage_thickn;
  G4double anode_to_SiPMs;
  
  G4double SiPM_between_long;
  G4double SiPM_between_short;
  G4int    SiPM_number;
  
  G4double Pb_box_lengt_xy;
  G4double Pb_box_lengt_z;
  G4double Pb_box_rel_pos;
  G4double Pb_box_z;
  

  //S1 AND S2 LENGTHS
  G4double drift_length;
  G4double el_length;
  
  G4double S2_rad;
  G4double S2_lenght;
  
  G4double S1_rad;
  G4double S1_lenght;
  
  //POSITIONS
  G4double vessel_z; 
  G4double anodeBracket_z;
  G4double anode_z;
  G4double gateBracket_z;
  G4double gate_z;
  G4double cathBracket_z;
  
  G4double teflon_cage_z;
  G4double TPB_tefloncage_z;
  G4double long_ring_z;
  G4double teflon_ring_z;
  G4double cath_ring_z;
  G4double ring0_z;
  
  G4double pmt_z;
  G4double plateUp_pmt_z;
  G4double enclosure_pmt_z;
  G4double plate_pmt_z;
  G4double PMTplateBottom1_pos_z;
  
  G4double drift_z;  
  G4double el_z;
  
  G4double S2_z;
  G4double S1_z;
  
  G4double encapsulation_z;
  
  G4double SiPMs_z;
  
};

field_cage_parameters version2_parameters();

G4LogicalVolume* get_world(field_cage_parameters const & fcp);
G4PVPlacement* GeometryV2();
G4PVPlacement* GeometryV2_TEST();

#endif
