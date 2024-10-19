'''
the finite difference method without the consideration of MSIS assumption
'''
import numpy as np 

class spherical_noCmax():
    def __init__(
        self,
        Faraday_constant = 96500, gas_constant = 8.314, temperature = 300,
        diffuse_solid = 1e-12, diffuse_liquid = 1e-6,
        concentration_liquid = 1e-1, 
        voltage_speed = 3e6, alpha = 0.5, k_zero = 1e-2, ele_num = 1, num = 200, r0 = 4e-6,
        voltage_init = 1, voltage_end = -1, voltage_standard = 0,
        deltaThetapoint = 2000, h = 1e-2, pomega_A =1.01 ,pomega_B =1.01, deltaX = 1e-4
    ):
        '''Faraday constant, gas constant, temperature'''
        self.Faraday_constant = Faraday_constant
        self.gas_constant = gas_constant
        self.temperature = temperature

        '''the diffusion coefficient of the electrolyte/electrode side'''
        self.diffuse_solid = diffuse_solid
        self.diffuse_liquid = diffuse_liquid

        '''the original concentration in electrolyte'''   
        self.concentration_liquid = concentration_liquid

        '''scan rate、alpha、k0, n, num, r0'''
        self.voltage_speed = voltage_speed                   
        self.alpha = alpha
        self.k = k_zero
        self.ele_num = ele_num
        self.num = num
        self.r0 = r0

        '''voltage'''
        self.voltage_init = voltage_init
        self.voltage_end = voltage_end
        self.voltage_standard = voltage_standard

        '''dimensionless'''
        self.theta_i = self.Faraday_constant/self.gas_constant/self.temperature*(self.voltage_init-self.voltage_standard)

        self.theta_v = self.Faraday_constant/self.gas_constant/self.temperature*(self.voltage_end-self.voltage_standard)

        self.dB = self.diffuse_solid/self.diffuse_liquid
        self.dA = self.diffuse_liquid/self.diffuse_liquid

        self.sigma = self.r0**2/self.diffuse_liquid*self.Faraday_constant/self.gas_constant/self.temperature*self.voltage_speed

        self.K0 = self.k*self.r0/self.diffuse_liquid

        '''grid parameter'''
        self.deltaThetapoint = deltaThetapoint
        self.deltaTheta = 2*abs(self.theta_i-self.theta_v) / deltaThetapoint
        self.h = h
        self.pomega_A = pomega_A
        self.pomega_B = pomega_B
        self.deltaX = deltaX
    
    def thomaxsolver(self,a,b,c,d):
        n = len(d)
        x = list(np.zeros(n))
        '''Calculate g_mod[]'''
        g_mod = [c[0] / b[0]]
        for j in range(1,n-1):
            gg = c[j] / (b[j] - g_mod[j-1]*a[j])
            g_mod.append(gg)
        '''Calculate d_mod[]'''
        d_mod = [d[0] / b[0]]
        for j in range(1,n):
            dd = (d[j] - d_mod[j-1]*a[j]) / (b[j] - g_mod[j-1]*a[j])
            d_mod.append(dd)
        '''Calculate x'''
        x[n-1] = d_mod[-1]  
        for j in range(n-2,-1,-1):
            x[j] = d_mod[j] - g_mod[j]*x[j+1]
        return x

    def grid(self):
        '''T grid'''
        self.deltaT = self.deltaTheta / self.sigma
        self.Tpoint = self.deltaThetapoint
        self.Tmax = 2*abs(self.theta_i-self.theta_v) / self.sigma

        '''X grid'''
        if self.dA > self.dB:
            dmax = self.dA
        else:
            dmax = self.dB
        Xmax = 6*np.sqrt(dmax*self.Tmax) + 1   # sphere diffuse

        self.XA = [1]
        deltaAh = self.h
        while self.XA[-1] < Xmax:
            self.XA.append(self.XA[-1]+deltaAh)
            deltaAh = deltaAh*self.pomega_A
        self.XApoint = int(len(self.XA))

        self.XB = [1]
        deltaBh = self.h
        while self.XB[-1] > 0:
            self.XB.append(self.XB[-1]-deltaBh)
            deltaBh = deltaBh*self.pomega_B
        self.XB = self.XB[:-1]
        self.XBpoint = int(len(self.XB))

        self.Xpoint = self.XApoint + self.XBpoint

    
    def fFunction(self,theta):
        a = np.exp(-1*self.alpha*theta)
        return a

    def calculate(self):
        self.grid()
        '''Thomas coefficients'''
        palpha_list = [-1,-1]
        pbeta_list = [-1,-1]
        pgamma_list = [-1,-1]

        for i in range(1,self.XApoint-1):
            '''Thomas coefficients for A'''
            RAhi = self.XA[i] - self.XA[i-1]
            RAhu = self.XA[i+1] - self.XA[i]

            palpha_list_a = -2*self.deltaT / (RAhi**2 + RAhi*RAhu) + 2 / self.XA[i] * (self.deltaT/(RAhi+RAhu))
            pbeta_list_a = 2*self.deltaT / (RAhu**2 +RAhi*RAhu) + 2*self.deltaT / (RAhi**2+RAhu*RAhi) + 1
            pgamma_list_a = -2*self.deltaT / (RAhu**2 + RAhi*RAhu) - 2 / self.XA[i] *(self.deltaT/(RAhi+RAhu))

            palpha_list.insert(0,pgamma_list_a)
            pbeta_list.insert(0,pbeta_list_a)
            pgamma_list.insert(0,palpha_list_a)

        for i in range(1,self.XBpoint-1):
            '''Thomas coefficients for B'''
            RBhi = self.XB[i] - self.XB[i-1]
            RBhu = self.XB[i+1] - self.XB[i]

            palpha_list_b = self.dB*(-2*self.deltaT / (RBhi**2 + RBhi*RBhu) + 2 / self.XB[i] * (self.deltaT/(RBhi+RBhu)))
            pbeta_list_b = self.dB*(2*self.deltaT / (RBhu**2 +RBhi*RBhu) + 2*self.deltaT / (RBhi**2+RBhu*RBhi)) + 1
            pgamma_list_b = self.dB*(-2*self.deltaT / (RBhu**2 + RBhi*RBhu) - 2 / self.XB[i] *(self.deltaT/(RBhi+RBhu)))

            palpha_list.append(palpha_list_b)
            pbeta_list.append(pbeta_list_b)
            pgamma_list.append(pgamma_list_b)

        
        palpha_list.insert(0,0)
        palpha_list.append(0)
        pbeta_list.insert(0,1)
        pbeta_list.append(1)
        pgamma_list.insert(0,0)
        pgamma_list.append(0)


        pdelta_list = list(np.zeros(self.XApoint+self.XBpoint))
        for i in range(self.XApoint):
            pdelta_list[i] = 1

        Theta = []
        J = []
        
        theta = self.theta_i
        for j in range(self.Tpoint):
            '''calculate value of theta at j'''
            if j < self.Tpoint/2 :
                theta = theta - self.deltaTheta
            else:
                theta = theta + self.deltaTheta
            Theta.append(theta)
            '''calculate thomas coefficient'''
            palpha_list[self.XApoint-1] = -1
            pbeta_list[self.XApoint-1] = 1 + self.h*self.fFunction(theta)*self.K0
            pgamma_list[self.XApoint-1] = -1*self.h*self.fFunction(theta)*self.K0*np.exp(theta) 
            pdelta_list[self.XApoint-1] = 0

            palpha_list[self.XApoint] = -1*self.h*self.fFunction(theta)*self.K0 / self.dB
            pbeta_list[self.XApoint] = 1 + self.h*self.fFunction(theta)*self.K0*np.exp(theta) / self.dB
            pgamma_list[self.XApoint] = -1
            pdelta_list[self.XApoint] = 0
            

            pdelta_list[self.Xpoint-1] = 0

            pdelta_list = self.thomaxsolver(palpha_list,pbeta_list,pgamma_list,pdelta_list)

            flux = -1*(-1*pdelta_list[self.XApoint-3]+4*pdelta_list[self.XApoint-2]-3*pdelta_list[self.XApoint-1]) / (2*self.deltaX)
            
            J.append(flux)
        return J,Theta

        
        
