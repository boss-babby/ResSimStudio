cd /home/user/ressim/research/exp
python3 e12d_dt.py /tmp/sd/cli_exp X:22:1:27.06:0.04438:50 2.5 all e12d_w25.json > e12d_w25.log 2>&1 &
python3 e12d_dt.py /tmp/sd/cli_exp X:22:0:32.68:0.03761:50 2.5 all e12d_n25.json > e12d_n25.log 2>&1
wait
python3 e12d_dt.py /tmp/sd/cli_exp X:22:1:27.06:0.04438:50 1.25 viol e12d_w125.json > e12d_w125.log 2>&1
echo E12D DONE
