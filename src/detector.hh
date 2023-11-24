#ifndef DETECTOR_H
#define DETECTOR_H

bool process_hits_anode(G4Step *step);
bool process_hits(G4Step *step);
bool process_hits_genratorCHECK(G4Step *step);

#endif
