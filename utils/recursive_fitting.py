import numpy as np
from concurrent.futures import ProcessPoolExecutor


def recursive_fitting(trj:np.ndarray, wlist:list, max_wokers:int=1):
    ### init fitting ###
    print('init fitting')
    superimpose = SuperImpose(trj, wlist, True, max_wokers)
    trj = superimpose()


    ### second fitting ###
    print('second fitting')
    superimpose = SuperImpose(trj, wlist, False, max_wokers)
    trj = superimpose()


    return trj



class SuperImpose():
    def __init__(self, trj:np.ndarray, wlist:list, whether_init:bool, max_wokers=1):
        self.trj = trj
        self.wlist = wlist
        self.n_frames = trj.shape[0]
        self.max_wokers = max_wokers

        if whether_init:
            self.reference_structure = trj[0]
        else:
            self.reference_structure = trj.mean(axis=0)
    

    def __call__(self):
        ### parallel process ###
        with ProcessPoolExecutor(max_workers=self.max_wokers) as executor:
            futures = []
            for i in range(self.n_frames):
                futures.append(executor.submit(self.rotate, i))

        return np.array([f.result() for f in futures])



    def rotate(self, target_index:int):
        target_structure = self.trj[target_index]

        ### matrix U ###
        U = np.empty((3,3))
        for i in range(3):
            for j in range(3):
                U[i][j] = sum( [self.wlist[n]*target_structure[n,i]*self.reference_structure[n,j] for n in range(target_structure.shape[0])] )

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
        
        self.print_progress(target_index)
        return fitted_structure


    def print_progress(self, i:int):
        if i % 10 == 0 or i+1 == self.n_frames:
            progress_frames = i+1
            progress_percentage = int(progress_frames/self.n_frames * 100)
            print(f"\rprogress: {progress_percentage}% ({progress_frames} frames)", end="")

        if i+1 == self.n_frames:
            print()