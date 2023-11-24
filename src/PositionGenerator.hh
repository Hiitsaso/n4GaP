#ifndef POSITIONGENERATOR_H
#define POSITIONGENERATOR_H

G4ThreeVector random_generator_inside_S1(std::optional<G4double> fixed_z);
G4ThreeVector random_generator_inside_S2(std::optional<G4double> fixed_x, std::optional<G4double> fixed_y, std::optional<G4double> fixed_z);
G4ThreeVector random_generator_tpb();

#endif
