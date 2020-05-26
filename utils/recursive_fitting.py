import numpy as np
from concurrent.futures import ProcessPoolExecutor
import time


def recursive_fitting(trj:np.ndarray, wlist:list, max_wokers:int, whether_recur:bool):

    ### init fitting ###
    print('First fitting into init structure')
    trj = do_fitting(trj, trj[0], wlist, max_wokers)

    ### second fitting ###
    if whether_recur:
        print('\nSecond fitting into mean structure')
        trj = do_fitting(trj, trj.mean(axis=0), wlist, max_wokers)

    print()
    return trj


def do_fitting(trj:np.ndarray, reference_structure:np.ndarray, wlist:list, max_wokers:int):
    myprocess = MyProcess(superimpose, trj.shape[0])
    with ProcessPoolExecutor(max_workers=max_wokers) as executor:
        futures = []
        for i in range(trj.shape[0]):
            futures.append(executor.submit(myprocess, i, trj[i], reference_structure, wlist))

    return np.array([f.result() for f in futures])



class MyProcess:
    def __init__(self, func, totalstep:int):
        self.totalstep = totalstep
        self.interval = totalstep // 100
        self.starttime = time.time()
        self.func = func
    
    def __call__(self, step:int, *args):
        ret = self.func(*args)

        if step%self.interval==0 or step==self.totalstep+1:
            progress = int((step+1) / self.totalstep * 100)
            elapsed_time = time.time() - self.starttime
            print(f'\rProgress: {progress}% {elapsed_time:.1f}s', end='')
                
        return ret




def superimpose(target_structure, reference_structure, wlist):
    ### matrix U ###
    U = np.empty((3,3))
    for i in range(3):
        for j in range(3):
            U[i][j] = sum( [wlist[n]*target_structure[n,i]*reference_structure[n,j] for n in range(target_structure.shape[0])] )

    ### matrix OMEGA ###
    OMEGA = np.empty((6,6))
    for i in range(3):
        for j in range(3):
            OMEGA[i][j] = 0
            OMEGA[i+3][j+3] = 0
            OMEGA[i+3][j] = U[i][j]
            OMEGA[i][j+3] = U.T[i][j]


    ### resolve Eigenvalue problem ###
    eig_val, eig_vec =np.linalg.eig(OMEGA)

    omegas = eig_vec.T

    ### split eig_vec to h and k ###
    hks = np.array([omegas[i] for i in range(6) if eig_val[i] > 0])
    hks = hks * np.sqrt(2) # root2

    k = np.empty((3,3))
    h = np.empty((3,3))
    for i in range(3):
        for j in range(3):
            k[i][j] = hks[j][i]
            h[i][j] = hks[j][i+3]
    

    ### rotation matrix R ###
    R = np.empty((3,3))
    for i in range(3):
        for j in range(3):
            R[i][j] = sum( [k[i][a]*h[j][a] for a in range(len(h))] )

    
    ### target_trj to fitted_trj ###
    fitted_structure = np.empty_like(target_structure)
    for i in range(3):
        for n in range(fitted_structure.shape[0]):
            fitted_structure[n,i] = sum( [R[i,j]*target_structure[n,j] for j in range(3)] )
    
    return fitted_structure
