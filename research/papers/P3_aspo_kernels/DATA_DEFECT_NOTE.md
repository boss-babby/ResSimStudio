# Draft note to the data authors (send before submitting P3)

**To:** corresponding author of Zang et al. (2024, ESSD 16:295–310) / GFZ Data Services (dataset 10.5880/GFZ.2.6.2023.004)
**Subject:** Possible timing defect in the `origin_time` column of the Äspö HRL AE catalogues

Dear colleagues,

Thank you for publishing the Äspö hydraulic-fracturing dataset (doi 10.5880/GFZ.2.6.2023.004, CC BY 4.0), which we are using in a
reanalysis of stage-level AE productivity.

While matching AE to the 1-s hydraulic records, we noticed that in both catalogue files we downloaded,
`seismic_data/AE_triggered/CATALOG` (196 events) and `seismic_data/AE_continuous/CATALOG` (19,756 events), the seconds field of
`origin_time` takes only the values 03, 04, 09, 11 and 12. These equal the day of month of each event. The `event_name` field, by
contrast, is monotonic and agrees with the hydraulic records to within the stage timing. Examples:

    continuous:  name 20150603093930.254880   origin_time 2015-06-03_09:39:03.254880
                 name 20150603152205.139680   origin_time 2015-06-03_15:22:03.139680
    triggered:   Origin_time 2015-06-03_10:04:03.352000   Name_in_conti_catalog 20150603100405.775680

In each case the minutes agree with the name, but the seconds equal the day (03).

Our analysis therefore derives the event times from `event_name`. Could you confirm whether this is a known formatting issue, and
whether a corrected version is planned? We will describe the issue neutrally in our manuscript and cite the dataset. We are happy to
share our parsing script.

Kind regards,
[authors]

---
*Internal: evidence in `research/exp/h6c_aspo.py` (parser), `h6c.log`, FINDINGS_LOG F39. Four names in the triggered/continuous
match file could not be parsed and were dropped (documented in P3 §2).*
