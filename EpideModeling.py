import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

st.title("Different Epidemic Models")

options = ["SIR", "SEIR", "SITR", "SLIAR", "SEQIJR"]
selectbox_selection = st.selectbox("Select Model", options)
st.write(f"Model selected is {selectbox_selection}")

if selectbox_selection == "SIR":
    
    st.latex(r'''\begin{align*}
             \frac{dS}{dt} &= -\beta S I \\
             \frac{dI}{dt} & = \beta S I -\alpha I \\
             \frac{dR}{dt} & = \alpha I
             \end{align*}''')
    
    alpha = st.number_input('Choose recovery rate')
    beta = st.number_input('Choose infection rate')
    
    # Simulate SIR Model
    def SIR():
    
        """time step in days ,infection rate, 
        recovery rate"""
        ndays = 100
        dt = .01         
        npts = int(ndays/dt)
    
        S = np.zeros(npts)  
        I = np.zeros(npts)   
        R = np.zeros(npts)    
        
        t = np.arange(npts)*dt
    
        #initial infective proportion
        I[0] = 0.001  
        
        #initial susceptible population
        S[0] = 1. - I[0]                           
        R[0] = 0.
    
        for i in range(npts-1):
            S[i+1] = S[i] - beta*(S[i]*I[i])*dt
            I[i+1] = I[i] + (beta*S[i]*I[i] - alpha*I[i])*dt
            R[i+1] = R[i] + (alpha*I[i])*dt
    
        fig = plt.figure(1); fig.clf()
        plt.plot(t, S,'r',lw=3,label='Susceptible')
        plt.plot(t, I,'g',lw=3,label='Infective')
        plt.plot(t, R,'b',lw=3,label='Removed')
    
        fig.legend(); plt.xlabel('Days'); plt.ylabel('Fraction of Population')
    
        Max_Infected = ((I).max()*100).round(2)
        Susceptible_Left = S[npts-1].round(4)*100
        data = [Max_Infected, Susceptible_Left]
        df = pd.DataFrame(data, columns=['Values'], index=['Max Infected',' Susceptible Left'])
    
        st.dataframe(df)
        st.pyplot(fig)
    SIR()

elif selectbox_selection == "SEIR":
    
    st.latex(r'''\begin{align*}
             \frac{dS}{dt} &= -\beta S I \\
             \frac{dI}{dt} & = \beta S I -\alpha I \\
             \frac{dR}{dt} & = \alpha I
             \end{align*}''')
    
    alpha = st.number_input('Choose recovery rate')
    beta = st.number_input('Choose infective rate')
    kappa = st.number_input('Choose rate at which they leave the exposed compartment')
    
    # Simulate SEIR Model
    def SEIR():
        ndays = 200          #number of days for the model
        dt = 0.01            #time step in days
        npts = int(ndays/dt) #number of data points (better resolution)
        
        #Making Placeholder Values
        S = np.zeros(npts)
        E = np.zeros(npts)
        I = np.zeros(npts)
        R = np.zeros(npts)
        
        t = np.arange(npts)*dt
        
        #Initial Infective Population
        I[0] = 0.0001
        
        #Initial Susceptible, Exposed, and Removed Population
        S[0] = 1. - I[0] 
        E[0] = 0.
        R[0] = 0.
        
        for i in range(npts-1):
            S[i+1] = S[i] - (beta*S[i]*I[i])*dt
            E[i+1] = E[i] + (beta*S[i]*I[i]-kappa*E[i])*dt
            I[i+1] = I[i] + (kappa*E[i]-alpha*I[i])*dt
            R[i+1] = R[i] + (alpha*I[i])*dt
            
        fig = plt.figure(1); fig.clf()
        plt.plot(t, S,'r',lw=2.5,label='Susceptible')
        plt.plot(t, E,'y',lw=2.5,label='Exposed')
        plt.plot(t, I,'g',lw=2.5,label='Infective')
        plt.plot(t, R,'b',lw=2.5,label='Removed')
    
        fig.legend(); plt.xlabel('Days'); plt.ylabel('Fraction of Population')
    
        Max_Infected = ((I).max()*100).round(2)
        Susceptible_Left = S[npts-1].round(4)*100
        data = [Max_Infected, Susceptible_Left]
        df = pd.DataFrame(data, columns=['Values'], index=['Max Infected',' Susceptible Left'])
    
        st.dataframe(df)
        st.pyplot(fig)
    SEIR()
    
elif selectbox_selection == "SITR":
    
    st.latex(r'''\begin{align*}
             \frac{dS}{dt} &= -\beta S I \\
             \frac{dI}{dt} & = \beta S I -\alpha I \\
             \frac{dR}{dt} & = \alpha I
             \end{align*}''')
    
    alpha = st.number_input('Choose recovery rate')
    beta = st.number_input('Choose infective rate')
    delta = st.number_input('Choose reduced infectivity rate')
    gamma = st.number_input('Choose rate where infectives get treated')
    eta = st.number_input('Choose rate of removal from treatment')

    
    # Simulate SITR Model
    def SITR():
        ndays = 250          #number of days for the model
        dt = 0.01            #time step in days
        npts = int(ndays/dt) #number of data points (better resolution)
    
        #Making Placeholder Values
        S = np.zeros(npts)
        I = np.zeros(npts)
        T = np.zeros(npts)
        R = np.zeros(npts)
        
        t = np.arange(npts)*dt
        
        #Initial Infective Population
        I[0] = 0.0001
        
        #Initial Susceptible, Treatment, and Removed Population
        S[0] = 1. - I[0] 
        T[0] = 0.
        R[0] = 0.
        
        for i in range(npts-1):
            S[i+1] = S[i] - (beta*S[i]*(I[i] + delta*T[i]))*dt
            I[i+1] = I[i] + (beta*S[i]*(I[i] + delta*T[i]) - (alpha + gamma)*I[i])*dt
            T[i+1] = T[i] + (gamma*I[i] - eta*T[i])*dt
            R[i+1] = R[i] + (alpha*I[i] + eta*T[i])*dt
            
        fig = plt.figure(1); fig.clf()
        plt.plot(t, S,'r',lw=2.5,label='Susceptible')
        plt.plot(t, I,'g',lw=2.5,label='Infective')
        plt.plot(t, T,'y',lw=2.5,label='Treated')
        plt.plot(t, R,'b',lw=2.5,label='Removed')
    
        fig.legend(); plt.xlabel('Days'); plt.ylabel('Fraction of Population')
        
        Max_Infected = ((I).max()*100).round(2)
        Susceptible_Left = S[npts-1].round(4)*100
        data = [Max_Infected, Susceptible_Left]
        df = pd.DataFrame(data, columns=['Values'], index=['Max Infected',' Susceptible Left'])
    
        st.dataframe(df)
        st.pyplot(fig)
    SITR()
    
elif selectbox_selection == "SLIAR":
    
    st.latex(r'''\begin{align*}
             \frac{dS}{dt} &= -\beta S I \\
             \frac{dI}{dt} & = \beta S I -\alpha I \\
             \frac{dR}{dt} & = \alpha I
             \end{align*}''')
    
    alpha = st.number_input('Choose recovery rate')
    beta = st.number_input('Choose infective rate')
    delta = st.number_input('Choose reduced infectivity rate')
    kappa = st.number_input('Choose rate where latent people either become infective or asymptomatic')
    rho = st.number_input('Choose fraction that goes to the infective')
    eta = st.number_input('Choose rate of removal from treatment')
    
    # Simulate SLIAR Model
    def SLIAR():
        ndays = 100          #number of days for the model
        dt = 0.01            #time step in days
        npts = int(ndays/dt) #number of data points (better resolution) 
        
        #Making Placeholder Values
        S = np.zeros(npts)
        L = np.zeros(npts)
        I = np.zeros(npts)
        A = np.zeros(npts)
        R = np.zeros(npts)
        
        t = np.arange(npts)*dt
        
        #Initial Infective Population
        I[0] = 0.0001
        
        #Initial Susceptible, Latent, Asymptomatic, and Removed Population
        S[0] = 1. - I[0] 
        L[0] = 0.
        A[0] = 0.
        R[0] = 0.
        
        for i in range(npts-1):
            S[i+1] = S[i] - (beta*S[i]*(I[i] + delta*A[i]))*dt
            L[i+1] = L[i] + (beta*S[i]*(I[i] + delta*A[i]) - kappa*L[i])*dt
            I[i+1] = I[i] + (rho*kappa*L[i] - alpha*I[i])*dt
            A[i+1] = A[i] + ((1-rho)*kappa*L[i] - eta*A[i])*dt
            R[i+1] = R[i] + (alpha*I[i] + eta*A[i])*dt
            
        fig = plt.figure(1); fig.clf()
        plt.plot(t, S,'r',lw=2.5,label='Susceptible')
        plt.plot(t, L,'c',lw=2.5,label='Latent')
        plt.plot(t, I,'g',lw=2.5,label='Infective')
        plt.plot(t, A,'y',lw=2.5,label='Asymptomatic')
        plt.plot(t, R,'b',lw=2.5,label='Removed')
    
        fig.legend(); plt.xlabel('Days'); plt.ylabel('Fraction of Population')
        
        Max_Infected = ((I).max()*100).round(2)
        Susceptible_Left = S[npts-1].round(4)*100
        data = [Max_Infected, Susceptible_Left]
        df = pd.DataFrame(data, columns=['Values'], index=['Max Infected',' Susceptible Left'])
    
        st.dataframe(df)
        st.pyplot(fig)
    SLIAR()
    
elif selectbox_selection == "SEQIJR":
    
    st.latex(r'''\begin{align*}
             \frac{dS}{dt} &= -\beta S I \\
             \frac{dI}{dt} & = \beta S I -\alpha I \\
             \frac{dR}{dt} & = \alpha I
             \end{align*}''')
    
    alphaI = st.number_input('Choose recovery rate from infective ')
    alphaJ = st.number_input('Choose recovery rate from isolation')
    beta = st.number_input('Choose infective rate')
    epsilonE = st.number_input('Choose reduced infectivity rate of exposed population')
    epsilonQ = st.number_input('Choose reduced contant rate of quarantine population')
    epsilonJ = st.number_input('Choose decreased infectivity rate of isolation population')
    gammaQ = st.number_input('Choose rate from exposed to quarantine')
    gammaJ = st.number_input('Choose rate from infective to isolation')
    kappaE = st.number_input('Choose rate from exposed to infective')
    kappaQ = st.number_input('Choose rate from quarantine to isolation')
    
    # Simulate SEQIJR Model
    def SEQIJR():
        ndays = 100          #number of days for the model
        dt = 0.01            #time step in days
        npts = int(ndays/dt) #number of data points (better resolution)
        
        #Making Placeholder Values
        S = np.zeros(npts)
        E = np.zeros(npts)
        Q = np.zeros(npts)
        I = np.zeros(npts)
        J = np.zeros(npts)
        R = np.zeros(npts)
        
        t = np.arange(npts)*dt
        
        #Initial Infective Population
        I[0] = 0.0001
        
        #Initial Susceptible, Exposed, Quarantine, Isolation, and Removed Population
        S[0] = 1. - I[0] 
        E[0] = 0.
        Q[0] = 0.
        J[0] = 0.
        R[0] = 0.
        
        for i in range(npts-1):
            S[i+1] = S[i] - (beta*S[i]*(epsilonE*E[i] + epsilonE*epsilonQ*Q[i] + I[i] + epsilonJ*J[i]))*dt
            E[i+1] = E[i] + (beta*S[i]*(epsilonE*E[i] + epsilonE*epsilonQ*Q[i] + I[i] + epsilonJ*J[i]) - (kappaE + gammaQ)*E[i])*dt
            Q[i+1] = Q[i] + (gammaQ*E[i] - kappaQ*Q[i])*dt
            I[i+1] = I[i] + (kappaE*E[i] - (alphaI + gammaJ)*I[i])*dt
            J[i+1] = J[i] + (kappaQ*Q[i] + gammaJ*I[i] - alphaJ*J[i])*dt
            R[i+1] = R[i] + (alphaI*I[i] + alphaJ*J[i])*dt
            
        fig = plt.figure(1); fig.clf()
        plt.plot(t, S,'r',lw=2.5,label='Susceptible')
        plt.plot(t, E,'y',lw=2.5,label='Exposed')
        plt.plot(t, Q,'c',lw=2.5,label='Quarantine')
        plt.plot(t, I,'g',lw=2.5,label='Infective')
        plt.plot(t, J,'m',lw=2.5,label='Isolation')
        plt.plot(t, R,'b',lw=2.5,label='Removed')
    
        fig.legend(); plt.xlabel('Days'); plt.ylabel('Fraction of Population')
        
        Max_Infected = ((I).max()*100).round(2)
        Susceptible_Left = S[npts-1].round(4)*100
        data = [Max_Infected, Susceptible_Left]
        df = pd.DataFrame(data, columns=['Values'], index=['Max Infected',' Susceptible Left'])
    
        st.dataframe(df)
        st.pyplot(fig)
    SEQIJR()