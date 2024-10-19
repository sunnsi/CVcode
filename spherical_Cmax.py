'''
the finite difference method under the MSIS assumption
'''
import numpy as np 

class spherical_Cmax():
    def __init__(
        self,
        Faraday_constant = 96485, gas_constant = 8.314, temperature = 300,
        diffuse_solid = 1e-12, diffuse_liquid = 1e-6,
        concentration_liquid = 1e-1, concentration_Cmax=1e0,
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
        self.concentration_Cmax = concentration_Cmax/concentration_liquid

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
        '''Thomas algorithm'''
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

            # liquid
            palpha_list_a = 2*self.deltaT*(1/(self.XA[i]*(RAhi+RAhu))-1/(RAhu*RAhi+RAhi*RAhi))
            pbeta_list_a = 1+2*self.deltaT*(1/(RAhu*RAhi+RAhu*RAhu)+1/(RAhu*RAhi+RAhi*RAhi))
            pgamma_list_a = -2*self.deltaT*(1/(RAhu*RAhi+RAhu*RAhu)+1/(self.XA[i]*(RAhi+RAhu)))

            palpha_list.append(palpha_list_a)
            pbeta_list.append(pbeta_list_a)
            pgamma_list.append(pgamma_list_a)


        for i in range(1,self.XBpoint-1):
            '''Thomas coefficients for B'''
            
            RBhi = self.XB[i] - self.XB[i-1]
            RBhu = self.XB[i+1] - self.XB[i]
            
            # solid
            palpha_list_b = 2*self.dB*self.deltaT*(1/(self.XB[i]*(RBhi+RBhu))-1/(RBhu*RBhi+RBhi*RBhi))
            pbeta_list_b = 1+2*self.dB*self.deltaT*(1/(RBhu*RBhi+RBhu*RBhu)+1/(RBhu*RBhi+RBhi*RBhi))
            pgamma_list_b = -2*self.dB*self.deltaT*(1/(RBhu*RBhi+RBhu*RBhu)+1/(self.XB[i]*(RBhi+RBhu)))

            palpha_list.insert(0,pgamma_list_b)
            pbeta_list.insert(0,pbeta_list_b)
            pgamma_list.insert(0,palpha_list_b)


        palpha_list.insert(0,0)
        palpha_list.append(0)
        pbeta_list.insert(0,1)
        pbeta_list.append(1)
        pgamma_list.insert(0,-1)
        pgamma_list.append(0)


        
        pdelta_list = list(np.zeros(self.XApoint+self.XBpoint))
        for i in range(self.XBpoint,self.XApoint+self.XBpoint):
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
            

            CA0 = pdelta_list[self.XBpoint-1]
            CB0 = pdelta_list[self.XBpoint]
            
            # electrode
            palpha_list[self.XBpoint-1] = 1
            
            pbeta_list[self.XBpoint-1] = -1-self.h/self.dB*self.fFunction(theta)*self.K0*(CB0/self.concentration_Cmax+np.exp(theta)) 

            pgamma_list[self.XBpoint-1] = self.h/self.dB*self.fFunction(theta)*self.K0*(1-CA0/self.concentration_Cmax) 
           
            # electrolyte
            palpha_list[self.XBpoint] = self.h*self.fFunction(theta)*self.K0*(CB0/self.concentration_Cmax+np.exp(theta))

            pbeta_list[self.XBpoint] = -1 + self.h*self.fFunction(theta)*self.K0*(CA0/self.concentration_Cmax-1)

            pgamma_list[self.XBpoint] = 1
            
            
            pdelta_list[self.XBpoint-1] = -self.h*self.fFunction(theta)*self.K0/self.dB/self.concentration_Cmax*CA0*CB0
            pdelta_list[self.XBpoint] = self.h*self.fFunction(theta)*self.K0/self.concentration_Cmax*CA0*CB0
            
            
            pdelta_list[0] = 0
            
            
            pdelta_list = self.thomaxsolver(palpha_list,pbeta_list,pgamma_list,pdelta_list)


            flux = -1*(-1*pdelta_list[self.XBpoint+2]+4*pdelta_list[self.XBpoint+1]-3*pdelta_list[self.XBpoint]) / (2*self.h)

            
            J.append(flux)

        return J,Theta

        
