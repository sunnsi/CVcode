'''fitness of Ipc data with the formula 
we deduced under the MSIS assumption'''
import numpy as np
import matplotlib.pyplot as plt
import xlrd
import xlwt

workbook=xlwt.Workbook(encoding='utf-8') 
idx=0

def funfit(r,nu,Ds,Cmax):
    '''
    Parameters
    ----------
    r : the radius of the spherical electrode
    nu : the scan rate
    Ds : the diffusion coefficient of the electrode side
    Cmax : MSIS
        
    Returns
    -------
    Ipc : the logarithm of peak current density
    '''
    temp = np.sqrt(Ds/nu)
    Ic = 1000*2.988*F*Cmax*r*nu*temp/(r+temp)
    Ipc = np.log10(Ic)
    return Ipc


book = xlrd.open_workbook('dataforMSIS.xls')
book1 = book.sheet_by_name('Sheet1')
speed = np.linspace(-4,6,11)
F = 96485
T = 300
Rg = 8.314


nu=np.power(10,speed)

plt.figure(figsize=(10,8),dpi=250)

shangC = ['red', 'blue', 'green']


'''
figure 1 of different Ds
data storage
'''
booksheet1=workbook.add_sheet('Sheet1', cell_overwrite_ok=True)
r = 1e-4
Cmax = [0.01,0.001,0.0001]
Ds =  [1e-7,1e-8,1e-9]
for k in range(len(Ds)):
    for i in range(len(Cmax)):
        Ipc3=[]
        for j in range(len(speed)):
            I3 = funfit(r,nu[j],Ds[k],Cmax[i])
            Ipc3.append(I3)
            booksheet1.write(j,idx,I3)
        idx=idx+1
        plt.plot(speed,Ipc3,color = shangC[i])
        
        
color1 = ['red', 'blue', 'green','red', 'blue', 'green','red', 'blue', 'green']
label1 = ['^','^','^','o','o','o','s','s','s']
for i in range(book1.ncols):
    a=book1.col_values(i)
    
    plt.scatter(speed,a,facecolor='none',edgecolor = color1[i],marker = label1[i],s=80)
  
plt.xlabel(r"$log(\nu)$",fontsize=30)
plt.ylabel(r"$log(I_{pc})$",fontsize=30)
plt.show()

'''
figue 2 of different r
data storage
'''
booksheet2=workbook.add_sheet('Sheet2', cell_overwrite_ok=True)
idx=0
book2 = book.sheet_by_name('Sheet2')
r = [1e-4,3e-4,10e-4]
plt.figure(figsize=(10,8),dpi=250)
for k in range(len(r)):
    for i in range(len(Cmax)):
        Ipc2=[]
        Ipc3=[]
        for j in range(len(speed)):
            I3 = funfit(r[k],nu[j],1e-7,Cmax[i])
            Ipc3.append(I3)
            booksheet2.write(j,idx,I3)
        idx=idx+1
        plt.plot(speed,Ipc3,color = shangC[i])
        
        
color2 = ['red','red','red', 'blue', 'blue','blue','green','green','green']
label2 = ['^','o','s','^','o','s','^','o','s']
for i in range(book2.ncols):
    a=book2.col_values(i)
    plt.scatter(speed,a,facecolor='none',edgecolor = color2[i],marker = label2[i],s=80)

plt.xlabel(r"$log(\nu)$",fontsize=30)
plt.ylabel(r"$log(I_{pc})$",fontsize=30)
plt.show()

workbook.save('datalineforMSIS.xls')