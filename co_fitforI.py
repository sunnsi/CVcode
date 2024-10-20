'''the differential evolution method to calculate the coefficient'''
import numpy as np
import matplotlib.pyplot as plt
import xlrd
import xlwt
from scipy.optimize import differential_evolution


def funfit(r,nu,Ds,Cmax):
    '''
    Parameters
    -------
    deduced formula in the article
    r : the radius of spherical electrode
    nu : the scan rate
    Ds : the diffusion efficient of the electrode side
    Cmax : MSIS
    
    Returns
    -------
    Ipc: the logarithm of peak current density
    '''
    temp2 = np.sqrt(Ds/nu)
    Ic = 1000*F*Cmax*nu*r*temp2/(r+temp2)
    return Ic

def fun1(theta,X,Y):
    y_pred = np.log10(theta*X)
    # residual sum of squares
    ss_res = np.sum((Y - y_pred) ** 2)
    # total sum of sqares
    ss_tot = np.sum((Y - np.mean(Y)) ** 2)
    r2 = 1-ss_res/ss_tot
    return r2

def fun2(theta,X,Y):
    return 1-fun1(theta[0],X,Y)

workbook=xlwt.Workbook(encoding='utf-8')
booksheet1=workbook.add_sheet('Sheet1', cell_overwrite_ok=True)

Ipc=np.array([])
x=np.array([])

'''data fetch'''
book = xlrd.open_workbook('dataforMSIS.xls')
book1 = book.sheet_by_name('Sheet1')
speed = np.linspace(-4,6,11)
F = 96485
T = 300
Rg = 8.314
nu=np.power(10,speed)

'''data of different diffusion coefficient Ds'''
r = 1e-4
Cmax = [0.01,0.001,0.0001]
Ds =  [1e-7,1e-8,1e-9]
for k in range(len(Ds)):
    for i in range(len(Cmax)):
        for j in range(len(speed)):
            para = funfit(r,nu[j],Ds[k],Cmax[i])
            x=np.append(x,para)

for i in range(book1.ncols):
    a=book1.col_values(i)
    a=np.array(a).reshape(-1)
    Ipc=np.append(Ipc,a)

'''data of different radius r'''
r = [1e-4,3e-4,10e-4]
book2 = book.sheet_by_name('Sheet2')

for k in range(len(Cmax)):
    for i in range(len(r)):
        for j in range(len(speed)):
            para = funfit(r[i],nu[j],1e-7,Cmax[k])
            x=np.append(x,para)

for i in range(book2.ncols):
    a=book2.col_values(i)
    a=np.array(a).reshape(-1)
    Ipc=np.append(Ipc,a)

Ipc10=np.power(10,Ipc)
xpc = np.log10(x)

'''differential_evolution method '''
bounds = [(0,10)]
result = differential_evolution(fun2,bounds,args=(x,Ipc),seed=1)
print("the optimal coeffcient is:",result.x[0])


gamma = result.x[0]
I_fit=np.log10(gamma*x)

'''data storage'''
i=0
for k in I_fit:
    booksheet1.write(i,0,k)
    i=i+1
i=0
for k in Ipc:
    booksheet1.write(i,1,k)
    i=i+1

test = np.linspace(-4,6,1000)
plt.figure(figsize=(10,8),dpi=250)
plt.plot(test,test,c='black')
plt.scatter(I_fit,Ipc,s=20,marker='o',c='w',edgecolors='g',label='$R^2=0.996$')
plt.xlabel(r"$log(\psi_{predict})$",fontsize=20)
plt.ylabel(r"$log(\psi_{actual})$",fontsize=20)

plt.legend(fontsize=20)
plt.show()


workbook.save('dataR2forI.xls')
