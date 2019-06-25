#!/bin/bash
rm -r /publicweb/m/mreid/iDM_AN_Plots/plots
cp -r plots /publicweb/m/mreid/iDM_AN_Plots/
python make_html_listing.py /publicweb/m/mreid/iDM_AN_Plots/plots/GenKinematics

find /publicweb/m/mreid/iDM_AN_Plots/plots -mindepth 0 -type d -exec cp plots/.htaccess {} \;


