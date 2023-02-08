import numpy as np
from random import seed
from random import sample

from f_observability import *


class Bus_system:
    def __init__(self, Ybus):
        self.Ybus=Ybus
        self.n_bus =  np.size(Ybus,1)
        self.n_branches = count_branches(Ybus)
        self.max_meas = self.n_branches*2 + self.n_bus

def count_branches(Ybus):
    n_branches = 0
    
    for line in Ybus:
        for i in line:
            n_branches += i
    n_branches = int(n_branches/2)
    return n_branches

def build_empty_measurement_plan(Ybus, max_med):
    meas_plan = np.zeros([max_med,7],np.int32)
    

    n_bus = len(Ybus)
    measurement = 0
    #medidas de fluxo
    for line in range(n_bus):
        for column in range(n_bus):
            if(Ybus[line][column] == 1):
                de = line+1
                para = column+1
                meas_plan[measurement][0] = measurement+1
                meas_plan[measurement][1] = de
                meas_plan[measurement][2] = para
                meas_plan[measurement][3] = 1
                meas_plan[measurement][4] = 1
                meas_plan[measurement][5] = de
                meas_plan[measurement][6] = 0
                
                measurement = measurement + 1

    #medidas de injeção
    for line in range(n_bus):
        bus = line+1
        meas_plan[measurement][0] = measurement + 1
        meas_plan[measurement][1] = bus
        meas_plan[measurement][2] = bus
        meas_plan[measurement][3] = 1
        meas_plan[measurement][4] = 2
        meas_plan[measurement][5] = bus
        meas_plan[measurement][6] = 0
        
        measurement = measurement + 1

    return meas_plan

def add_random_measurement(meas_plan, possible_meas):
    chosen_meas = sample(possible_meas, 1)
    possible_meas = [bus for bus in possible_meas if bus not in chosen_meas]

    if meas_plan[chosen_meas[0]-1,6] == 0 : 
        meas_plan[chosen_meas[0]-1,6] = 1 
        return 1
    else :
        return 0

def remove_desactivated_measurements(meas_plan): #arquivo apenas com as medidas ativadas
    medidas_desativadas = []
    num_medidas_desativadas = 0
    cont = 0
    for line in meas_plan:
        if line[6] == 0:
            medidas_desativadas.append(cont)
            num_medidas_desativadas += 1
        cont = cont + 1

    #Deleta medidas da lista
    meas_plan = np.delete(meas_plan,medidas_desativadas, axis=0)
    #Corrige contagem inicial
    for i in range(meas_plan.shape[0]):
        meas_plan[i,0] = i + 1
    return meas_plan

def get_terminal_bus(Ybus):
    (n_bus, _) = Ybus.shape
    terminal_buses = []
    for bus in range(n_bus):
        if (Ybus[bus,:] == 1).sum(axis=0) == 1:
            terminal_buses.append(bus+1)
    
    return terminal_buses

def get_initial_meas_plan(meas_plan, terminal_buses):
    for bus in terminal_buses:
        for line in range(meas_plan.shape[0]):
            if (meas_plan[line,1] == bus+1):
                meas_plan[line,6] = 1

    return meas_plan

def get_most_relevant_buses(Ybus):
    n_neighbors = []
    i = 1
    for bus in Ybus:
        n_neighbors.append((i,sum(bus)))
        i= i + 1
    
    n_neighbors.sort(key=lambda x: x[-1], reverse = True) 
    relevant_buses = [i[0] for i in n_neighbors] 
    n_relevants = int(len(Ybus)/4)

    return(relevant_buses[0:n_relevants])

if __name__ == '__main__':
    Ybus_File = "Rede_Polonia.txt"
    Ybus = np.loadtxt(Ybus_File, dtype='i', delimiter=',')

    get_most_relevant_buses(Ybus)
    # terminal_buses = get_terminal_bus(Ybus)

    # power_system = Bus_system(Ybus)

    # empty_plan = build_empty_measurement_plan(Ybus, power_system.max_meas)
    # pre_processed_meas = get_initial_meas_plan(empty_plan, terminal_buses)
    print("Acabou")