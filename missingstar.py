import matplotlib.pyplot as plt
import numpy as np
from numpy import random
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.time import Time
from datetime import datetime, timedelta
from astropy.coordinates import SkyCoord
import pandas as pd
from tqdm import tqdm, trange

#read_data는 list 주소를 받아 별들의 좌표와 HR 번호에 대응되는 별 이름을 출력한다.
def read_data(path):
    coor=[]
    data=pd.read_csv(path,delimiter="\t")
    data=np.array(data)
    star_name=[' ']*10000
    for i in tqdm(data,desc='read the data and convert it into array'):
        l=len(i[3].split(' '))
        if l==3:
            number,name,_,mag,spec=int(i[0]),i[1],i[2],float(i[5]),i[6]
            ra_h,ra_m,ra_s=map(float,i[3].split(' '))
            dec_d,dec_m,dec_s=map(float,i[4][1:].split(' '))
            ra=ra_h*15.+ra_m*0.25+ra_s/240.
            dec=dec_d+dec_m/60.+dec_s/3600.
            if i[4][0]=='-':
                dec*=-1
            coor.append([ra,dec,mag,number])
            star_name[number]=name
    return np.array(coor),np.array(star_name)

#sidereal 함수는 location, UT time을 받으면 sidereal time을 출력해준다. 
def sidereal(lat,lon,time,timezone):
    ut=timedelta(hours=timezone)
    obs_site=EarthLocation(lat=lat*u.deg,lon=lon*u.deg)
    obs_time=Time(time-ut,scale='utc',location=obs_site)
    LST=obs_time.sidereal_time('mean')
    return LST.hour

# 각 별들의 sinAcosa, cosAcosa, sina값을 갖고 있다. limmag1, limmag2는 각각 그래프에 표시되는 한계등급과 missing star 등급
def riseset(location,data,limmag1,limmag2,missing_num,add_num):
    lat,lon,time,timezone=location[0],location[1],location[2],location[3]
    LST=sidereal(lat,lon,time,timezone)
    data=data.T #data는 행이 순서대로 ra,dec,mag,numb가 됨.
    ra,dec,mag,numb=data[0],data[1],data[2],data[3]
    LST_deg=LST*15.
    hour_angle=LST_deg-ra

    c=np.pi/180.
    
    row1=np.sin(hour_angle*c)*np.cos(dec*c)
    row2=np.cos(hour_angle*c)*np.cos(dec*c)
    row3=np.sin(dec*c)

    A=np.array([row1,row2,row3,mag,numb])
    transf_mat=np.array([[1.,0.,0.,0.,0.],[0,np.sin(lat*c),-np.cos(lat*c),0.,0.],[0,np.cos(lat*c),np.sin(lat*c),0.,0.],[0.,0.,0.,1,0.],[0.,0.,0.,0.,1]])
    
    B=transf_mat@A

    alt=(np.arcsin(1)-np.arcsin(B[2]))/np.arcsin(1)
    sina=B[2]
    cosa=np.cos(np.arcsin(B[2]))
    new_row1=alt*(B[0])/cosa
    new_row2=alt*(B[1])/cosa

    B[0],B[1]=new_row1,new_row2

    C1=(B[2]>0) #고도가 0보다 큰가?
    C2=(mag<limmag2) #random missing star를 빼기 위한 등급조건 
    C3=np.logical_and(mag<limmag1,mag>=limmag2) #sky chart에 표시할 별의 등급조건



    B1=B[:,np.logical_and(C1,C2)].T #random missing star 조건을 만족하는 별
    B2=B[:,np.logical_and(C1,C3)].T #위 조건에 해당하지 않는 별

    lis=list(range(len(B1)))
    a=random.choice(lis,size=missing_num,replace=False)
    nota=np.array([x for x in lis if x not in a])
    
    missing_star=B1[a]
    B1=B1[nota]
    B=np.concatenate([B1,B2])

    #additional star 만들기
    phi,theta=random.rand(add_num)*(2.*np.pi),random.rand(add_num)*(np.pi/2)
    rand_mag=random.rand(add_num)*3.5-1.   #-1.0mag~2.5mag
    X,Y,Z=np.cos(phi)*np.cos(theta),np.sin(phi)*np.cos(theta),np.sin(theta)
    dum_num=np.array([-1]*add_num)
    add_star=np.array([X,Y,Z,rand_mag,dum_num])


    return missing_star.T,np.array(B).T,add_star

#random으로 성도를 회전하는 역할
def rotate(star_list,theta):
    R=np.array([[np.cos(theta),np.sin(theta),0,0,0],[-np.sin(theta),np.cos(theta),0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])
    new_starlist=R@star_list
    return new_starlist


def graphplot(missing_star,star_list,add_star,name,pro_num):
    draw_circle1=plt.Circle((0,0),1.,fill=False)
    draw_circle2=plt.Circle((0,0),1.,fill=False)
    #solution part
    solpath='./Problem/Sol_'+str(pro_num).zfill(3)+'.png' #정답지 파일 주소
    fig1,ax1=plt.subplots(1,1,figsize=(20,20))
    ax1.scatter(missing_star[0],-missing_star[1],10**(2.5-0.4*missing_star[3]),c='r')
    ax1.scatter(star_list[0],-star_list[1],10**(2.5-0.4*star_list[3]),c='k')
    ax1.scatter(add_star[0],-add_star[1],10**(2.5-0.4*add_star[3]),c='b')

    nametag=missing_star[4].astype(np.int64)
    for i in range(len(nametag)):
        ax1.text(missing_star[0][i]-.05,-missing_star[1][i]-.05,name[nametag[i]])
    ax1.set_title('Solution #'+str(pro_num),fontsize=15)
    ax1.add_artist(draw_circle1)
    ax1.axis('off')
    fig1.savefig(solpath,bbox_inches='tight',transparent=False)
    #problem part
    propath='./Problem/Pro_'+str(pro_num).zfill(3)+'.png' #문제지 파일 주소
    fig2,ax2=plt.subplots(1,1,figsize=(20,20))
    ax2.scatter(star_list[0],-star_list[1],10**(2.5-0.4*star_list[3]),c='k')
    ax2.scatter(add_star[0],-add_star[1],10**(2.5-0.4*add_star[3]),c='k')
    ax2.set_title('Problem #'+str(pro_num),fontsize=15)
    ax2.add_artist(draw_circle2)
    ax2.axis('off')
    fig2.savefig(propath,bbox_inches='tight',transparent=False)


def main():
    star_path='./bsc.tsv'  #별 목록의 주소
    data,name=read_data(star_path)

    #원하는 사진 장수, missing star, additional star 개수
    inp,missing_num,add_num=100,3,3
    

    for i in trange(inp,desc='Generate the test paper'):
        rand_hour=random.randint(24)
        rand_minute=random.randint(60)
        rand_time=datetime(2022,10,12,rand_hour,rand_minute,00)
        location=[37.5,127.05,rand_time,9]
        missing_star,star_list,add_star=riseset(location,data,5.5,2.5,3,3)
        theta=(random.rand(1)*(2*np.pi))[0]
        missing_star,star_list,add_star=rotate(missing_star,theta),rotate(star_list,theta),rotate(add_star,theta)
        graphplot(missing_star,star_list,add_star,name,i) #problem, solution 주소는 이 함수 내에서 바꾸기.
    
    print("It's done!")

main()