# TIGERISS_Abundances
2026 TIGERISS Abundance Updates

This project is to make updated predictions of the measurements of the Trans-Iron Galactic Element Recorder for the International Space Station (TIGERISS) in the ultra-heavy galactic cosmic ray (UHGCR) regime. We combine a prediction model developed for the CALorimetric Electron Telescope (CALET) with simulations run using GEometry ANd Tracking 4 (Geant4). Predictions are made for 1 year of operation aboard the International Space Station.

## Previous Model
Current TIGERISS \cite{Rauch:2023apa} observation predictions build off of an older CALET prediction model \cite{RAUCH20141444} which uses differential geometry factors as a function of East-West inclination angle to account for geomagnetic screening against energy spectra derived from other cosmic ray experiments. Instrument geometry factors are calculated by reading in the dimensions and positions of all detector layers. Starting from the smallest geometric layer, random points are chosen and trajectories are drawn in both directions. Events that hit all layers are binned by incidence angle ($\theta$) and hit orientation relative to East-West angle. This returns both a 1-D differential geometry factor dependent on $\theta$ and a 2-D differential geometry factor incorporating orientation relative to the East-West angle, which is used in the Stoermer Approximation \cite{Rauch:2023apa}.

Top-of-instrument (TOI) elemental counts are determined using integrated spectra derived from ACE-CRIS measurements at solar maximum and minimum and HEAO-3-C2 measurements at higher energies. These events are then filtered  using ray-traced geomagnetic cutoff rigidities at ISS orbital altitude and weighted by the time spent across all geomagnetic coordinates. The TIGERISS trigger is modeled by an SSD energy deposit equivalent to  half of that deposited by a vertically incident <sub>5</sub>B particle. 

The TOI threshold energies for events to trigger in all layers as a function of charge and angle are determined through species-specific Geant4 simulations and used as an additional constraint on detection. The higher of either the local geomagnetic rigidity cutoff or the instrument triggering threshold establishes the lower integration bound used to evaluate final elemental abundance predictions. TIGERISS analysis will enforce charge consistency between triggers, so Nilsen cross sections \cite{Nilsenetal:1995} are used to determine survival fraction for all events that are improperly counted.


## Completed Updates

### Triggering Thresholds
Triggering threshold energies have been generated using Geant4 which simulates millions of events passing through the instrument and records the energy (MeV/nuc) deposited in each layer. The energy at which a majority of events trigger the detector is set as the TOI threshold. The thresholds are binned by angle and incident particle charge for Z=5-82 and $\theta$=0-70. Updated energy thresholds were created for the TIGERISS 9x9 SSD model placed aboard the ISS Starboard Overhead X-Direction (SOX) module with mature geometry and material compositions. These results were presented at the 2025 International Cosmic Ray Conference in poster "Geant4 Simulations of Geometry Factor and Interaction and Energy Losses for TIGERISS".

### Updated Material Interactions
TIGERISS is proposed to match SuperTIGER's measurements within one year aboard the ISS. Spallation within the instrument will occur at rates described by Nilsen et al. Early prediction models examined the effects of SuperTIGER materials on "survival fraction", the fraction of events not undergoing a charge-changing interaction, for all charges that TIGERISS will measure using analytic formulas. This fraction is incorporated into the abundance predictions.

Updating the prediction model to find the survival fractions of particles passing through the mature TIGERISS instrument will show the effects of transitioning from balloon flight to space travel and the use of SSDs. Material densities and chemical compositions determined by the TIGERISS team are used to calculate the updated survival fraction for the mature instrument.

### Accounting for Geomagnetic Orientation
The ISS will fly between -51.6 °N and 51.6 °N as it orbits the Earth, passing all longitudes in a 90-minute orbit. As the ISS changes geographic latitude, TIGERISS will rotate about its central axis, resulting in orientation changes. Additionally, each coordinate has its own magnetic declination, defined as the difference between True North (the North Pole) and the magnetic North Pole.  Combining ISS orientations relative to True North, which are a function of latitude, and magnetic declinations, which are determined for each latitude and longitude using the 2025 World Magnetic Model (WMM) at a constant altitude of 450 km, yields the instrument's geomagnetic heading. 

TIGERISS differential geometry factors are calculated for rotations, allowing us to see the effect of the major axis being unaligned with geomagnetic north. The fraction of time spent at each geomagnetic heading is determined using ISS orbital data. The most time is spent at the highest and lowest altitudes. We transition from using two geometry factors, one with the major axis aligned with the geomagnetic North Pole and one with the minor axis aligned with the geomagnetic North Pole, to using time-weighted differential geometry factors for all geomagnetic headings.


## Updates in Progress

### Field-of-View Obstructions

### Implementing a New Rigidity Model

### Updating Input Spectra

### Accounting for Space Weather Effects

