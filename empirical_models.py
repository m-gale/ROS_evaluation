# -*- coding: utf-8 -*-
"""
Author: Matthew Gale matthew.gale@anu.edu.au
Empirical ROS model calculations for ROS model evaluation

Inputs:
rh - relative humidity (%)
temp - temperature (C)
sf_hazard - surface fuel hazard
FL_s - surface fuel load (tha)
nsf_hazard - near-surface fuel hazard
nsf_height - near-surface fuel height (cm)
FL_ns - near-surface fuel load (tha)
ef_height - elevated fuel height (m)
ef_hazard - elevated fuel hazard 
tfl - total fuel load (tha)
ws_net - 10 m open windspeed (ms-1)
waf - wind adjustment factor
df - drought factor
ffdi - forest fire danger index


"""

#%%

def calculate_vesta1_afdrs(rh, temp, nsf_hazard, sf_hazard, ws_net, nsf_height):

    #moisture content
    mc=2.76+0.124*rh-0.0187*temp
    #moisture function
    mf=18.35*mc**-1.495

    ws_kph=ws_net*60*60/1000
    
    if ws_kph>5:
        ros_vesta1=(30+1.531*((ws_kph-5)**0.858)*(sf_hazard**0.93)*((nsf_hazard*nsf_height*100)**0.637)*1.03)*mf
    else:
        ros_vesta1=30*mf
            
    return ros_vesta1


#%%

def calculate_vesta2(waf, df, FL_s, FL_ns, ef_height, ef_hazard, rh, temp, ws_net):
        
    #wind speed adjustment
    ws_kph=ws_net*60*60/1000
    u2=ws_kph/waf
    
    #moisture content
    mc=2.76+0.124*rh-0.0187*temp
    if mc<4.1:
        md=1
    else:
        if mc>24:
            md=0
        else:
            md=0.9082 + (0.1206*mc) - (0.03106*mc**2) + (0.001853*mc**3) - (0.00003467*mc**4)
    fa=1.008/(1 + 104.9 * math.exp(-0.9306*df))
    M=md*fa
    
    #phase 1 transitions
    if (FL_s+FL_ns)<1:
        p2=0
    else:
        gx_2 = -23.9315 + 1.7033 *u2 + 12.0822 * M + 0.95236 * (FL_s+FL_ns)
        p2=1/(1+math.exp(-gx_2))
    
    #phase 2 transitions
    #requires undestorey fuel height
    hu= -0.1 + 0.06 * ef_hazard + 0.48 * ef_height
    #requires ros2
    r2=0.19591*(u2**0.8257) * (((FL_s+FL_ns)/10)**0.4672) * hu**0.495 * M
    gx_3 = -32.3074 + 0.2951 * ws_kph + 26.8734 * M
    if r2<0.3:
        p3=0
    else:
        p3=1/(1+math.exp(-gx_3))
        
    #r1 abd r3
    if u2 > 2:
        r1=0.03+0.05024 * ((u2-1)**0.92628) * (((FL_s+FL_ns)/10)**0.79928) * M
    else:
        r1=0.03*M           
    r3=0.05235 * (ws_kph**1.19128) * M
    
    #overall ROS
    if p2<0.5:
        ros_vesta2 = r1 * (1-p2) + r2 * p2
    else:
        ros_vesta2= r1 * (1-p2) + r2 * p2 * (1-p3) +r3*p3

    return ros_vesta2

#%%

def calculate_mk5(ffdi, tfl):

    
    ros_mk5=0.0012*ffdi*tfl

            
    return ros_mk5

#%%

def calculate_10wind(ws_net):

    
    ros_10wind=0.1*(ws_net*60*60)

            
    return ros_10wind

#%%