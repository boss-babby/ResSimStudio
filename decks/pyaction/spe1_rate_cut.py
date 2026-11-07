# SPE1 PYACTION example (ResSim C37 / OPM Flow embedded Python)
# The script runs at every call; module globals persist between calls.
import opm_embedded

summary_state = opm_embedded.current_summary_state
schedule = opm_embedded.current_schedule

if 'setup_done' not in locals():
    setup_done = True
    done = False
    fopt_limit = 1.0e7          # stb

if not done and summary_state['FOPT'] > fopt_limit:
    done = True
    days = summary_state.elapsed() / 86400.0
    opm_embedded.OpmLog.info('FOPT %.4g stb at day %.1f: producer rate cut to 5000 stb/d (WBHP %.1f psia)'
                             % (summary_state['FOPT'], days, summary_state.well_var('PROD', 'WBHP')))
    schedule.insert_keywords("""
WCONPROD
 'PROD' OPEN ORAT 5000 1* 1* 1* 1* 1000 /
/
""")
