'''fitness of Ipc data with traditional formulas derived from Bard's book'''
import numpy as np
import matplotlib.pyplot as plt
import xlrd
import xlwt

workbook=xlwt.Workbook(encoding='utf-8') 
idx=0

def funbook(r,nu,Cmax):
    '''
    Parameters
    -------
    the traditional formula for peak current density Ipc in Bard's book
    r : the radius of spherical electrode
    nu : the scan rate
    Cmax : MSIS
    
    Returns
    -------
    Ipc: the logarithm of peak current density
    '''
    Dl = 1e-5
    Clstar = 1e-3
    temp = 0.399*np.sqrt(F*F*F/Rg/T*Dl*nu)*Clstar+F*Dl*Clstar/r*0.912
    Ipc = np.log10(temp)+3
    return Ipc

book = xlrd.open_workbook('datafornoMSIS.xls')
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
        Ipc2=[]
        for j in range(len(speed)):
            I2 = funbook(r,nu[j],Cmax[i])
            Ipc2.append(I2)
            booksheet1.write(j,idx,I2)
        idx=idx+1
        plt.plot(speed,Ipc2,color = shangC[i])
        
        
color1 = ['red', 'blue', 'green','red', 'blue', 'green','red', 'blue', 'green','red', 'blue', 'green']
label1 = ['^','^','^','^','o','o','o','o','s','s','s','s']
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
r = [1e-4,3e-4,10e-4]
book2 = book.sheet_by_name('Sheet2')
plt.figure(figsize=(10,8),dpi=250)
for k in range(len(r)):
    for i in range(len(Cmax)):
        Ipc2=[]
        for j in range(len(speed)):
            I2 = funbook(r[k],nu[j],Cmax[i])
            Ipc2.append(I2)
            booksheet2.write(j,idx,I2)
        idx=idx+1
        plt.plot(speed,Ipc2,color = shangC[k])
        
color2 = ['red', 'blue','green','red','red','red', 'blue', 'blue','blue','green','green','green']
label2 = ['^','o','s','^','o','s','^','o','s','^','o','s']
for i in range(book2.ncols):
    a=book2.col_values(i)
    plt.scatter(speed,a,facecolor='none',edgecolor = color2[i],marker = label2[i],s=80)
plt.xlabel(r"$log(\nu)$",fontsize=30)
plt.ylabel(r"$log(I_{pc})$",fontsize=30)
plt.show()

workbook.save('datalinefornoMSIS.xls')