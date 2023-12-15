#include <n4-all.hh>

//#include <G4SystemOfUnits.hh>

#include "detector.hh"
#include "GeometryV2.hh"

bool process_hits_anode(G4Step *step){			
	const std::string& filename_map = "Test_anode_vessel_2.txt";
	//~ const std::string& filename_map_primaries = "nan.txt";
	G4Track* track = step -> GetTrack();
	G4String particleType = track->GetDefinition()->GetParticleName();			
	G4int hits_check = 0;
	
	G4ThreeVector position     = track -> GetVertexPosition();
	G4ThreeVector position_hit = step  -> GetPostStepPoint() -> GetPosition();
	int colWidth = 20;
	
	//G4cout << "*************************************  :)  " << particleType << "  (:  *************************************"  << G4endl;   
	
	if (particleType == "opticalphoton"){ 
	track -> SetTrackStatus(fStopAndKill);
	
	hits_check = 1;
	//G4cout << "*************************************  :)  OUCH  (:  *************************************"  << G4endl; 
				
	}
	
	std::ofstream file1(filename_map, std::ios::app);
	file1 << std::left << std::setw(colWidth) << position.x();
	file1 << std::left << std::setw(colWidth) << position.y();
	file1 << std::left << std::setw(colWidth) << position.z();
	file1 << std::left << std::setw(colWidth) << position_hit.x();
	file1 << std::left << std::setw(colWidth) << position_hit.y();
	file1 << std::left << std::setw(colWidth) << position_hit.z();
	file1 << std::left << std::setw(colWidth) << hits_check;
	file1 << std::endl;
	file1.close();
	
	return true;
}

bool process_hits(G4Step *step){	
    field_cage_parameters fcp = version2_parameters();
				
	//~ const std::string& filename_map = "files_4SiPM/OP_Ar_128nm_S2_Z5mm_0mm_5000000.txt";
	const std::string& filename_map = "/home/hiitsaso/n4_GaProyect/nan.txt";
	//~ const std::string& filename_map = "/home/hiitsaso/n4_GaProyect/files/SiPMs/yes_WLS_no_EL/from_S2/OP_Ar_128nm_S2_10000000_fixed_S2z_4.txt";
	//~ const std::string& filename_map = "/home/hiitsaso/n4_GaProyect/files/SiPMs/yes_WLS_no_EL/from_S2/OP_Ar_128nm_S2_10000000_random_v3.txt";
	//~ const std::string& filename_map = "/home/hiitsaso/n4_GaProyect/files/SiPMs/yes_WLS_no_EL/from_S2/test/OP_Ar_128nm_S2_10000000_random_getEnergy_withTPBextra_withoutOpticalSurfaceTPBtoVessel.txt";
	//~ const std::string& filename_map = " files/no_WLS_no_EL/from_S2/OP_Ar_128nm_S2_10000000_random.txt";
	//~ const std::string& filename_map = "files_OneSiPM/v3_MESH_RAD/test.txt";
	//~ const std::string& filename_map = "files_OneSiPM/v3_MESH_RAD/random_test6.txt";
	G4Track* track = step -> GetTrack();
	G4String particleType = track->GetDefinition()->GetParticleName();			
	G4int hits_check = 0;
	
	G4ThreeVector position     = track -> GetVertexPosition();
	G4ThreeVector position_hit = step  -> GetPreStepPoint() -> GetPosition();
	G4double energy            = step  -> GetPreStepPoint() -> GetTotalEnergy();
	
	int colWidth = 20;
	
	auto zmin    = fcp.ring0_z - fcp.ring0_length/2; 
    auto zmax    = fcp.ring0_z + fcp.ring0_length/2; 
    auto xy_lim  = fcp.ring1_rad;
    
	//~ if (position.x() >= xy_lim && position.x() <= -xy_lim  &&  position.y() >= xy_lim && position.y() <= -xy_lim  && position.z() >= zmax && position.z() <= zmin){ 
	//no se como conseguir que no tenga en cuenta las que se generan ene sa parte random.
    	
		if (particleType == "opticalphoton"){ 
		track -> SetTrackStatus(fStopAndKill);
	
		hits_check = 1;
		//~ G4cout << "*************************************  :)  OUCH  (:  *************************************"  << G4endl; 
				
		}
	
		std::ofstream file1(filename_map, std::ios::app);
		file1 << std::left << std::setw(colWidth) << position.x(); 
		file1 << std::left << std::setw(colWidth) << position.y();
		file1 << std::left << std::setw(colWidth) << position.z();
		file1 << std::left << std::setw(colWidth) << position_hit.x();
		file1 << std::left << std::setw(colWidth) << position_hit.y();
		file1 << std::left << std::setw(colWidth) << position_hit.z();
		file1 << std::left << std::setw(colWidth) << energy;
		file1 << std::left << std::setw(colWidth) << hits_check;
		file1 << std::endl;
		file1.close();
		
		//~ if(file1.is_open()){ 
		//~ G4cout << "*************************************  :)  FILE DONE  (:  *************************************"  << G4endl; 
		//~ }else{
		//~ G4cerr << "Error opening file: " << strerror(errno) << G4endl;
		//~ G4cerr << "Error: File does not exist - " << filename_map << G4endl;
		//~ G4cout << "*************************************  :)  FILE NOOOOT DONE  (:  *************************************"  << G4endl; } 
				
	//~ }
	
	return true;
}


bool process_hits_genratorCHECK(G4Step *step){	
    field_cage_parameters fcp = version2_parameters();
				
	const std::string& filename_map = "files_OneSiPM/OpticalPhoton_420nm_S1_5mm_1000000_generatorCHECK.txt";
	G4Track* track = step -> GetTrack();
	G4String particleType = track->GetDefinition()->GetParticleName();			
	G4int hits_check = 0;
	
	G4ThreeVector position     = track -> GetVertexPosition();
	G4ThreeVector position_hit = step  -> GetPreStepPoint() -> GetPosition();
	
	int colWidth = 20;
	
	auto zmin    = fcp.ring0_z - fcp.ring0_length/2; 
    auto zmax    = fcp.ring0_z + fcp.ring0_length/2; 
    auto xy_lim  = fcp.ring1_rad;
    	
		if (position == position_hit){ 
		//~ track -> SetTrackStatus(fStopAndKill);
	
		hits_check = 1;
		//~ G4cout << "*************************************  :)  OUCH  (:  *************************************"  << G4endl; 
				
		}
	
		std::ofstream file1(filename_map, std::ios::app);
		file1 << std::left << std::setw(colWidth) << position.x(); 
		file1 << std::left << std::setw(colWidth) << position.y();
		file1 << std::left << std::setw(colWidth) << position.z();
		file1 << std::left << std::setw(colWidth) << position_hit.x();
		file1 << std::left << std::setw(colWidth) << position_hit.y();
		file1 << std::left << std::setw(colWidth) << position_hit.z();
		file1 << std::left << std::setw(colWidth) << hits_check;
		file1 << std::endl;
		file1.close();
	
	return true;
}
