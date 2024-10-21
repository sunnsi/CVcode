'''data generation without the consideration of MOS constraint'''
import spherical_noCmax as pn
import spherical_Cmax as pc

import numpy as np
import matplotlib.pyplot as plt
import xlwt

workbook=xlwt.Workbook(encoding='utf-8')

def plotCVtotal1():
    '''
    Spehrical electrode with MOS constraint and without MOS constraint.
    Different diffusion coefficient Ds.
    '''
    x = np.linspace(-4,6,11)
    voltage_speed = np.power(10,x)
    
    Ds_Value = [1e-7,1e-8,1e-9]
    shangDs = ['$10^{-7}$','$10^{-8}$','$10^{-9}$']
    shangC = ['purple','red', 'blue', 'green']
    shangS = ['^', '8',  's']
    plt.figure(figsize=(10,8),dpi=250)
    
    '''data storage'''
    booksheet1=workbook.add_sheet('Sheet1', cell_overwrite_ok=True)
    m=0

    for j in range(len(Ds_Value)):
        I_max = []
        for i in range(len(voltage_speed)):
            cv = pn.spherical_noCmax(
                    Faraday_constant = 96485, gas_constant = 8.314, temperature = 300,
                    diffuse_solid = Ds_Value[j], diffuse_liquid = 1e-5,
                    concentration_liquid = 1e-3,
                    voltage_speed =voltage_speed[i], alpha = 0.5, k_zero = 1e0, ele_num = 1, num = 200, r0 = 1e-4,
                    voltage_init = 0.5, voltage_end = -0.5, voltage_standard = 0, deltaThetapoint = 200, h = 1e-4, 
                    pomega_A =1.1 ,pomega_B =1.01, deltaX = 1e-4)
        
            J, theta = cv.calculate()
            
            mylabel = r'$D_s$={a} $cm^2/s $'.format(a=shangDs[j])
            mymarker = r"{a}".format(a=shangS[j])
            Imax=np.max(-np.array(J))
            # Normalize dimensionless peak current density
            Imax = 1000*96485*cv.concentration_liquid*cv.diffuse_liquid/cv.r0*Imax
            phi = np.log10(Imax)
            booksheet1.write(i,m,phi)
            I_max.append(phi)
        m+=1
        plt.plot(x,I_max, marker=mymarker, color=shangC[0],label = mylabel)

    plt.xlabel(r"$log(\nu)$",fontsize=25)
    plt.ylabel(r"$log(I_{pc})$",fontsize=25)
    plt.legend(bbox_to_anchor=(1, 1.15),ncol=4)
    plt.show()
    
plotCVtotal1()


def plotCVtotal2():
    '''
    Spehrical electrode with MOS constraint and without MOS constraint,
    Different radius R.
    '''
    x = np.linspace(-4,6,11)
    voltage_speed = np.power(10,x)
    
    R_Value = [1e-4,3e-4,10e-4]
    shangDs = ['$1*10^{-4}$','$3*10^{-4}$','$1*10^{-3}$']
    shangC = ['purple','red', 'blue', 'green']
    shangS = ['^', '8',  's']
    plt.figure(figsize=(10,8),dpi=250)
    
    '''data storage'''
    booksheet2=workbook.add_sheet('Sheet2', cell_overwrite_ok=True)
    m=0
    for j in range(len(R_Value)):
        I_max = []
        for i in range(len(voltage_speed)):
            cv = pn.spherical_noCmax(
                    Faraday_constant = 96485, gas_constant = 8.314, temperature = 300,
                    diffuse_solid = 1e-7, diffuse_liquid = 1e-5,
                    concentration_liquid = 1e-3,
                    voltage_speed =voltage_speed[i], alpha = 0.5, k_zero = 1e0, ele_num = 1, num = 200, r0 = R_Value[j],
                    voltage_init = 0.5, voltage_end = -0.5, voltage_standard = 0, deltaThetapoint = 200, h = 1e-4, 
                    pomega_A =1.1 ,pomega_B =1.01, deltaX = 1e-4)
        
            J, theta = cv.calculate()
            
            
            Imax=np.max(-np.array(J))
            # Normalize dimensionless peak current density
            Imax = 1000*96485*cv.concentration_liquid*cv.diffuse_liquid/cv.r0*Imax

            mylabel = r'$R_0$={a} $cm$'.format(a=shangDs[j])
            mymarker = r"{a}".format(a=shangS[j])
            # print(Imax)
            phi = np.log10(Imax)
            booksheet2.write(i,m,phi)
            I_max.append(phi)
        m+=1
        plt.plot(x,I_max, marker=mymarker, color=shangC[0],label = mylabel)

    plt.xlabel(r"$log(\nu)$",fontsize=25)
    plt.ylabel(r"$log(I_{pc})$",fontsize=25)
    plt.legend(bbox_to_anchor=(1, 1.15),ncol=4)
    plt.show()
    
plotCVtotal2()
workbook.save('datafornoMOS.xls')