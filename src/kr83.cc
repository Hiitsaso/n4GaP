#include <n4-all.hh>

#include <G4SystemOfUnits.hh>
#include <G4RandomDirection.hh>

#include "G4IonTable.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "kr83.hh"
#include "GeometryV2.hh" //creo que se puede quitar
#include "PositionGenerator.hh"


void kr83_generator(G4Event* event, G4double energy_32_ , G4double energy_9_ , G4double probGamma_9_, G4double lifetime_9_) {
    //auto gamma = nain4::find_particle("gamma");
                                      
    // From the TORI /ENSDF data tables.
    

    auto energy_Xrays_      = n4::scale_by(keV , {1.383, 1.435, 1.580, 1.581, 1.632, 1.647 , 1.699, 1.707, 1.906, 1.907, 12.405  , 12.598, 12.651, 14.104, 14.111, 14.231  , 14.311, 14.326 });
    auto probability_Xrays_ = n4::scale_by(0.01, {0.099, 0.060, 0.21 , 1.9  , 1.1  , 0.0094, 0.09 , 0.14 , 0.008, 0.025, 3.90E-05,  5.05 ,  9.8  ,  0.70 ,  1.36 ,  0.00429,  0.179,  0.0064});
    // Cumulative sum
    for (size_t i=1; i < probability_Xrays_.size(); i++) { probability_Xrays_[i] += probability_Xrays_[i-1];}

    // Set particle type searching in particle table by name
    //auto particle_defgamma_ = G4ParticleTable::GetParticleTable()->
    //  FindParticle("gamma");
    //auto particle_defelectron_ = G4ParticleTable::GetParticleTable()->
    //  FindParticle("e-");
      
    auto particle_defgamma_ = nain4::find_particle("gamma");
    auto particle_defelectron_ = nain4::find_particle("e-");
      
    // Ask the geometry to generate a position for the particle
    G4ThreeVector position = random_generator_inside_S1({});

   // First transition (32 kEv) Always one electron. Set it's kinetic energy.
   // Decide if we emit an X-ray..
   const double probXRay = G4UniformRand();
   double eKin32 = energy_32_;
   double eXray = 0.;
   if (probXRay < probability_Xrays_[probability_Xrays_.size() - 1]) {
    size_t kSel = probability_Xrays_.size() - 1;
    for (size_t k = 1; k != probability_Xrays_.size(); k++) {
      if ((probXRay >= probability_Xrays_[k - 1]) && (probXRay < probability_Xrays_[k])) {
        kSel = k - 1;
        break;
      }
    }
    eKin32 = energy_32_ - energy_Xrays_[kSel];
    eXray = energy_Xrays_[kSel];
  }

  G4double mass = particle_defelectron_->GetPDGMass();
  G4double energy = eKin32 + mass;
  G4double pmod = std::sqrt(energy * energy - mass * mass);
  G4ThreeVector momentum_dir32 = G4RandomDirection();
  G4double px = pmod * momentum_dir32.x();
  G4double py = pmod * momentum_dir32.y();
  G4double pz = pmod * momentum_dir32.z();
  // Set starting time of generation
  G4double time = 0;
  
  G4PrimaryVertex* vertex = new G4PrimaryVertex(position, time);
  
  // create new primaries and set them to the vertex
  
  G4PrimaryParticle* particle1 = new G4PrimaryParticle(particle_defelectron_);
  particle1->SetMomentum(px, py, pz);
  particle1->SetPolarization(0., 0., 0.);
  particle1->SetProperTime(time);
  vertex->SetPrimary(particle1);
  
  G4PrimaryParticle* particle2;

  if (eXray > (0.0001 * keV)) {
    pmod = eXray;
    G4ThreeVector momentum_dirXr = G4RandomDirection();
    px = pmod * momentum_dirXr.x();
    py = pmod * momentum_dirXr.y();
    pz = pmod * momentum_dirXr.z();
    particle2 = new G4PrimaryParticle(particle_defgamma_);
    particle2->SetMomentum(px, py, pz);
    particle2->SetPolarization(0., 0., 0.);
    particle2->SetProperTime(time);
    vertex->SetPrimary(particle2);
  } 
    // set the second gamma. Decide first if we have a 9.4 keV gamma,
    // on an EC
    
    const double probGam9 = G4UniformRand();
    G4ThreeVector momentum_dir9 = G4RandomDirection();
    //
    // finite lifetime for the 9 keV line
    //
    const double time9 = lifetime_9_ * G4RandExponential::shoot();
    
    G4PrimaryParticle* particle3;

    if (probGam9 < probGamma_9_) { // a soft gamma (i.e., an X-ray.. )
        energy = energy_9_ + mass;
        pmod = std::sqrt(energy * energy - mass * mass);
        G4ThreeVector momentum_dir9 = G4RandomDirection();
        px = pmod * momentum_dir9.x();
        py = pmod * momentum_dir9.y();
        pz = pmod * momentum_dir9.z();
        particle3 = new G4PrimaryParticle(particle_defelectron_);
        particle3->SetMomentum(px, py, pz);
        particle3->SetProperTime(lifetime_9_ * G4RandExponential::shoot());
        vertex->SetPrimary(particle3);
  } else { // a soft electron. (Electron Conversion )
     // Calculate cartesian components of momentum for the most energetic EC
      energy = energy_9_ + mass;
      pmod = std::sqrt(energy*energy - mass*mass);
      px = pmod * momentum_dir9.x();
      py = pmod * momentum_dir9.y();
      pz = pmod * momentum_dir9.z();
      G4PrimaryParticle* particle3 =
      new G4PrimaryParticle(particle_defelectron_);
      particle3->SetMomentum(px, py, pz);
      particle3->SetProperTime(time9);
      vertex->SetPrimary(particle3);
   }
   
    event->AddPrimaryVertex(vertex);
    
    //Delete particle events (not working, segmentation fault)
    
}
