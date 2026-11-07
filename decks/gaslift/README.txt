Gas-lift optimisation decks (C35, ECLIPSE / OPM LIFTOPT, WLIFTOPT, GLIFTOPT)
Source: OPM/opm-tests gaslift/ (Open Database License, ODbL). Unchanged except GASLIFT-04: PATHS 'MODEL05' -> 'include'.

  GASLIFT-01.DATA   9x9x2, 6 THP producers (THP 60, LRAT 5E3, WEFAC 0.8) + 4 VREP water injectors, VFPPROD 1 with an
                    ALQ axis 0..4E5 sm3/d. Groups PLAT-1 (PLAT-1A, PLAT-1B), PLAT-2; FIELD ORAT 10E3.
                    Schedule: LIFTOPT 4E4 (increment), all OP* optimised; 2021-08-01..2021-10-01 lift gas off (WLIFTOPT
                    NO 0.0); 2022-01-01 FIELD NONE + minimum economic gradient 6E-3; 2022-04-01 weight 3 for OP-A01 /
                    OP-B01 / OP-C02; 2022-08-01 FIELD ORAT 10E3 again + GLIFTOPT PLAT-1 maximum lift gas 9E5.
                    The wells do not flow without lift gas (FOPT 0.54E6 sm3 without, 9.70E6 with optimisation).
  GASLIFT-04.DATA   three-phase MODEL05 PVT (include/), VFPPROD with ALQ type GRAT (include/well_vfp.ecl), METRIC.

Reference: the commercial simulator results plotted in opm-tests gaslift/plots/ECL (no summary files are published;
values digitised from the plots, about +-1 %). OPM Flow has no reference run for these decks.

                      ResSim engine     commercial simulator
  GASLIFT-01  FOPT    9.70E6 sm3        ~10.25E6  (-5.3 %; late-time water distribution differs, wells die earlier)
              FWPT    0.74E6            ~1.09E6
              FPR     120.4 bar         ~111
  GASLIFT-04  FOPT    1.2517E7 sm3      ~1.255E7  (-0.3 %)
              FGPT    4.415E8 sm3       ~4.42E8   (-0.1 %)
              FWPT    4.20E5            ~4.5E5
              FPR     119.0 bar         ~116.5

Tests: ressim_tests -glift (import, lift vs no lift, increments, WLIFTOPT NO, GLIFTOPT group limit, economic gradient,
restart, round trip, GASLIFT-04 vs the commercial simulator).
Method: ENGINES_AND_MATH 6.30; user guide section 7j.
