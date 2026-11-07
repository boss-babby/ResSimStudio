IMPORTANT!

1) The sampling rates differ depending on the experiment

Sampling rate HF1-HF3: 0.2s
Sampling rate HF4-HF6: 0.01s

2) The timing precision is 1sec for HFs and 1 min for the Lungeon tests
	depending on the precision of the starting time of the record.
	--> The relative timing is correct but the starting time could deviate.
	--> Record starting times where reconstructed for HF4 to HF6 to our best knowledge

3) The hydraulic data of experiment HF6 is split into two files beccause it was conducted on two days.
The refracturing stages were stimulated one day later then the initial fracturing stage.

FILE FORMAT
(semicolon seperated)

HF1 -HF3:
unix_time;year-month-day;hour-min-sec.sec;flow_rate;interval_pressure;packer_pressure

HF4 -HF6:
unix_time;year-month-day;hour-min-sec.sec;flow_rate;interval_pressure;p1,p3

unix_time 			in seconds since 1970-01-01 00:00:00.000
flow rate			in L/min
interval pressure 	in MPa
packer pressure 	in MPa
p1,p3				in MPa

