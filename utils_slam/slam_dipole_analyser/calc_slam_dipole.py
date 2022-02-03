#!/bin/python

import sys,math
    
class slam_dipole_mod:

    def __init__(self,*args):   
        ''' 
            args ...

            pos_integral / config / LP_MO / type_info / MM_CNT / QM_CNT / Opti_Done
        
        '''

        self.pos_integral = float(args[0])                      # position integral val
        self.config = args[1]                                   # main configuration
        self.mo     = args[2]                                   # molecular orbital info
        self.type_info  = args[3]                               #.type_info info

        self.mm_cnt = int(args[4])                              # MM_CNT
        self.qm_cnt = int(args[5])                              # QM_CNT

        self.if_success = int(args[6])              # if SLAM optimisation is fully done

        self.mm_config = [[ 0 for i in range(6)] for j in range(self.mm_cnt)]   # 1,2,3,4,5,6 - name.type_info charge x y z
        self.qm_config = [[ 0 for i in range(6)] for j in range(self.qm_cnt)]   # 1,2,3,4     - name core_q sp_q x y z
        self.mo_config = [[ 0 for i in range(6)] for j in range(self.qm_cnt)]   # 1,2,3,4     - s, px, py, pz  (2) (3) (4) (5) tags

        # VECTOR WRITE PARAMETERS

        self.v_mod=1
        self.v_rho=4        # ARROW SCALING ?

        self.v_radius = 0.24
        self.v_rgb_r = 0
        self.v_rgb_g = 0
        self.v_rgb_b = 255      # 0 0 255 -> solid blue
        self.v_tail_flag = "1"

        self.sumofatomicdip=1       # WHAT DOES THIS PART DO?
    
        try:
            with open(self.config,"r") as c:            # READ CONFIG INFO

                rl = c.readline()

                for i in range(self.mm_cnt):            # READ MM CONFIG
                                        # LIST ROW STRIDE 6 ... NAME / TYPE(c/s) / CHARGE / x / y / z
                    rl = c.readline()
                    spl= rl.split()

                    self.mm_config[i][0] = spl[0]
                    self.mm_config[i][1] = spl[1]
                    #self.mm_config[i][2] = 
                    self.mm_config[i][3] = float(spl[2])
                    self.mm_config[i][4] = float(spl[3])
                    self.mm_config[i][5] = float(spl[4])

                for i in range(self.qm_cnt):

                    rl = c.readline()
                    spl= rl.split()
                    self.qm_config[i][0] = spl[0]
                    self.qm_config[i][3] = float(spl[1])
                    self.qm_config[i][4] = float(spl[2])
                    self.qm_config[i][5] = float(spl[3])

        except FileNotFoundError:
            print(FileNotFoundError)

        try:                            # READ TYPE INFO ... SAVE CHARGE INFO
            with open(self.type_info,"r") as t:

                rl = t.readline()
        
                for i in range(self.mm_cnt):

                    rl = t.readline()
                    spl= rl.split()

                    if spl[0] == self.mm_config[i][0] and spl[1] == self.mm_config[i][1]:
                        self.mm_config[i][2] = float(spl[2])

                for i in range(self.qm_cnt):

                    rl = t.readline()
                    spl= rl.split()

                    if spl[0] == self.qm_config[i][0]:

                        self.qm_config[i][1] = float(spl[1])
                        self.qm_config[i][2] = float(spl[2])

        except FileNotFoundError:
            print(FileNotFoundError)

        try:                            # READ MO INFO ... 
            with open(self.mo,"r") as m:
                                    # LIST ROW STRIDE 6 ... NAME / ENERGY(eV) / s / px / py / pz
                rl = m.readline()

                for i in range(self.qm_cnt):

                    rl = m.readline()
                    spl= rl.split()

                    self.mo_config[i][0] = spl[0]           # mo atom name
                    self.mo_config[i][1] = float(spl[1])    # mo atom energy
                    self.mo_config[i][2] = float(spl[2])    # mo atom s_
                    self.mo_config[i][3] = float(spl[3])    # mo atom px
                    self.mo_config[i][4] = float(spl[4])    # mo atom py
                    self.mo_config[i][5] = float(spl[5])    # mo atom pz

        except FileNotFoundError:
            print(FileNotFoundError)

    def get_base(self):

        if self.if_success == 1:
            self.base_if_optimised = True
        else:
            self.base_if_optimised = False

        self.base_pos_int = self.pos_integral
        self.base_if_charge_neutral = None
        
        self.base_mm_charge_sum = 0.
        self.base_qm_charge_sum = 0.
        self.base_total_charge_sum = 0.

        for i in range(self.mm_cnt):
            self.base_mm_charge_sum += self.mm_config[i][2]
        for i in range(self.qm_cnt):
            self.base_qm_charge_sum += (self.qm_config[i][1] + self.qm_config[i][2])

        self.base_total_charge_sum = self.base_mm_charge_sum + self.base_qm_charge_sum
    
        # GET IF SYSTEM CAHRGE NEUTRAL
        if self.base_total_charge_sum == 0.:
            self.base_if_charge_neutral = True
        else:
            self.base_if_charge_neutral = False

        #Input Convention
        #
        #self.mm_config = [[ 0 for i in range(6)] for j in range(self.mm_cnt)]   # 1,2,3,4,5,6 - name.type_info charge x y z
        #self.qm_config = [[ 0 for i in range(6)] for j in range(self.qm_cnt)]   # 1,2,3,4     - name core_q sp_q x y z
        #self.mo_config = [[ 0 for i in range(6)] for j in range(self.qm_cnt)]   # 1,2,3,4     - s, px, py, pz 

    def get_mm_dipole(self):

        mm_dip = [[ 0. for i in range(3) ] for j in range(self.mm_cnt)]

        scale=1.    # SCALE TEST FOR THE SIZE OF VIS DIPOLES IN *.vesta File .. default 1.

        for i in range(self.mm_cnt):    # GET MM DIPOLE ... VALUES MAY CHANGE DEPENDING ON THE GIVEN REFERENCE FRAME ( IF CHARGE IS NOT NEUTRAL )
        # simply ...
        # u = r*q

            mm_dip[i][0] = scale*self.mm_config[i][3]*self.mm_config[i][2]          # get dip x elem
            mm_dip[i][1] = scale*self.mm_config[i][4]*self.mm_config[i][2]          # get dip y elem
            mm_dip[i][2] = scale*self.mm_config[i][5]*self.mm_config[i][2]          # get dip z elem

        return mm_dip

    def get_qm_dipole(self):        # GET QM DIPOLE ... VALUES MAY CHANGE DEPENDING ON THE GIVEN REFERENCE FRAME ( IF CHARGE IS NOT NEUTRAL )

        qm_dip = [[ 0. for i in range(3) ] for j in range(self.qm_cnt)]

        scale=1.    # SCALE TEST FOR THE SIZE OF VIS DIPOLES IN *.vest File .. default 1.
    
        for i in range(self.qm_cnt):
        #   dip_elem =  r*(qc + qs) + <psi| r | psi>*qs
        # u = rc*qc + <psi| rs |psi>*qs
        # u = rc*qc + <psi| rcs + rc |psi>*qs
        # u = rc*qc + rc*qs + <psi| rcs |psi>*qs

            qm_dip[i][0]=scale*self.qm_config[i][3]*(self.qm_config[i][1] + self.qm_config[i][2]) \
                            + 2.*self.mo_config[i][2]*self.mo_config[i][3]*self.pos_integral*self.qm_config[i][2]
            qm_dip[i][1]=scale*self.qm_config[i][4]*(self.qm_config[i][1] + self.qm_config[i][2]) \
                            + 2.*self.mo_config[i][2]*self.mo_config[i][4]*self.pos_integral*self.qm_config[i][2]
            qm_dip[i][2]=scale*self.qm_config[i][5]*(self.qm_config[i][1] + self.qm_config[i][2]) \
                            + 2.*self.mo_config[i][2]*self.mo_config[i][5]*self.pos_integral*self.qm_config[i][2]

        return qm_dip

    def get_cluster_dipole(self):   # SUM OF THE CLUSTER DIPOLE MOMENT

        self.mm_dip = self.get_mm_dipole()
        self.qm_dip = self.get_qm_dipole()
    
        self.cluster_dip = [ 0. for i in range(3) ]     # Total Cluster Dipole

        self.cluster_dip_qm = [0. for i in range(3) ]   # Part QM   (Density)
        self.cluster_dip_mm = [0. for i in range(3) ]   # Part MM   (PointCharge)

        for i in range(self.mm_cnt):
            self.cluster_dip[0] += self.mm_dip[i][0]
            self.cluster_dip[1] += self.mm_dip[i][1]
            self.cluster_dip[2] += self.mm_dip[i][2]
            
            self.cluster_dip_mm[0] += self.mm_dip[i][0]
            self.cluster_dip_mm[1] += self.mm_dip[i][1]
            self.cluster_dip_mm[2] += self.mm_dip[i][2]

        for i in range(self.qm_cnt):
            self.cluster_dip[0] += self.qm_dip[i][0]
            self.cluster_dip[1] += self.qm_dip[i][1]
            self.cluster_dip[2] += self.qm_dip[i][2]

            self.cluster_dip_qm[0] += self.qm_dip[i][0]
            self.cluster_dip_qm[1] += self.qm_dip[i][1]
            self.cluster_dip_qm[2] += self.qm_dip[i][2]

    def write_intro(self):

        print("#"*90)
        print("")
        print(" Sp-Lone pair involved Atomistic Model ( S L A M ) ")
        print("")
        print(" (1) Toolkit - Dipole Moment Calculator")
        print("")
        print(" Compatibility   : SLAM Version Above 2.2 ")
        print(" (Last edited    : 07 - 11 - 2021 W.K. Jee)")
        print("")
        print("#"*90)
        print("")


    def write(self):            # WRITE GENERAL OUTPUT

    # write '.vesta' dipole vector file

        with open("dip_slam.vesta","w") as f:

            f.write("#VESTA_FORMAT_VERSION 3.3.0\n")
            f.write("MOLECULE\n")
            f.write("TITLE\n")
            f.write(" SLAM_DIPOLE\n")
            f.write("STRUC\n")      # STRUCTURE

            cnt=0
            offset_cnt=0
            for i in range(self.mm_cnt):

                if self.mm_config[i][1] == "c":
                    cnt += 1
                    offset_cnt += 1
                    atom_label = self.mm_config[i][0] + str(cnt)
                    f.write("%3d%6.3s%10.4s%12s" % (cnt,self.mm_config[i][0],atom_label,"1.0000"))
                    f.write("%12.6f%12.6f%12.6f" % (self.mm_config[i][3],self.mm_config[i][4],self.mm_config[i][5]))
                    f.write("%12s%12.2s\n" % ("1","-"))
                    f.write("%21.6f%12.6f%12.6f%12.2f\n" % (0.,0.,0.,0.))

            cnt=0
            for i in range(self.qm_cnt):
                cnt += 1
                atom_label = self.qm_config[i][0] + str(cnt)
                f.write("%3d%6.3s%10.4s%12s" % (cnt+offset_cnt,self.qm_config[i][0],atom_label,"1.0000"))
                f.write("%12.6f%12.6f%12.6f" % (self.qm_config[i][3],self.qm_config[i][4],self.qm_config[i][5]))
                f.write("%12s%12.2s\n" % ("1","-"))
                f.write("%21.6f%12.6f%12.6f%12.2f\n" % (0.,0.,0.,0.))
        
            f.write("    0 0 0 0 0 0 0\n")  # END FLAG
            # END OF WRITING STRUC      

            # BOND INFO
            if self.mm_cnt != 0:
                f.write("SBOND\n")
                f.write("%3d%6.3s%6.3s%70s\n" % (1,self.mm_config[0][0],self.qm_config[0][0],"0.00000    2.82146  0  1  1  0  1  0.250  2.000 127 127 127"))
                f.write("    0 0 0 0\n")

            # WRITE VECTOR
            f.write("VECTR\n")

            v_cnt=0
            for i in range(self.mm_cnt):

                if self.mm_config[i][1] == "c":
                    v_cnt += 1  # if core i.e., no vis atomic dipole ... pass
                if self.mm_config[i][1] == "s":

                    # r_shell - r_core = r_cs
                    rx = self.mm_config[i][3] - self.mm_config[i-1][3]
                    ry = self.mm_config[i][4] - self.mm_config[i-1][4]
                    rz = self.mm_config[i][5] - self.mm_config[i-1][5]

                    ux = rx*self.mm_config[i][2]*self.v_mod     # ~.mm_config[x][2] -> charge
                    uy = ry*self.mm_config[i][2]*self.v_mod
                    uz = rz*self.mm_config[i][2]*self.v_mod

                    ux=ux*-1.
                    uy=uy*-1.
                    uz=uz*-1.               # To make the arrow toward shell direction (only if shellcharge < 0.)
                    '''
                    amp = math.sqrt(ux*ux + uy*uy + uz*uz)
                    scl = math.exp(-amp/self.v_rho)             # ????
                    ux = ux*scl
                    uy = uy*scl
                    uz = uz*scl
                    '''
                    f.write("%3d%15.6f%12.6f%12.6f%12.1d\n" % (v_cnt,ux,uy,uz,0))
                    f.write("%3d" % (v_cnt))
                    f.write("    0 0 0 0\n")
                    f.write(" 0 0 0 0 0\n")

            for i in range(self.qm_cnt):

                v_cnt += 1
                ux = 2.*self.mo_config[i][2]*self.mo_config[i][3]*self.pos_integral*self.qm_config[i][2]*self.v_mod
                uy = 2.*self.mo_config[i][2]*self.mo_config[i][4]*self.pos_integral*self.qm_config[i][2]*self.v_mod
                uz = 2.*self.mo_config[i][2]*self.mo_config[i][5]*self.pos_integral*self.qm_config[i][2]*self.v_mod

                ux=ux*-1.
                uy=uy*-1.
                uz=uz*-1.               # To make the arrow toward shell direction (only if shellcharge < 0.)
                '''
                amp = math.sqrt(ux*ux + uy*uy + uz*uz)
                scl = math.exp(-amp/self.v_rho)
                ux = ux*scl
                uy = uy*scl
                uz = uz*scl
                '''
                f.write("%3d%15.6f%12.6f%12.6f%12.1d\n" % (v_cnt,ux,uy,uz,0))
                f.write("%3d" % (v_cnt))
                f.write("    0 0 0 0\n")
                f.write(" 0 0 0 0 0\n")
            f.write("    0 0 0 0 0\n")
            f.write("    0 0 0 0 0\n")
            # END WRITING VECTOR

            f.write("VECTT\n")
            # VECTOR FORMAT
            for i in range(v_cnt):
                f.write("%3d%9.3f%12d%12d%12d%12.3s\n" % (i+1,self.v_radius,self.v_rgb_r,self.v_rgb_g,self.v_rgb_b,self.v_tail_flag))
            f.write(" 0 0 0 0 0\n")

        self.write_intro()
        print(" Configuration / Atomic Dipole Moment Elements (Defined at chosen atom's core position)")
        print("")
        print("-"*90)
        print(" Species.  x           y           z               u_x           u_y           u_z")
        print("-"*90)


        self.atomic_dip_sum = [0. for i in range(3)]

        for i in range(self.mm_cnt):

            if   self.mm_config[i][1] == "c":

                atom_name   = self.mm_config[i][0]
                pos_x       = self.mm_config[i][3]
                pos_y       = self.mm_config[i][4]
                pos_z       = self.mm_config[i][5]              # if chosen atom is 'shel' then get 'core' position

                print("%3s%2s%10.6f%12.6f%12.6f" % (atom_name,"c",pos_x,pos_y,pos_z))

            elif self.mm_config[i][1] == "s":

                if self.sumofatomicdip == 1:
                    atom_name   = self.mm_config[i][0]
                    pos_x       = self.mm_config[i][3]
                    pos_y       = self.mm_config[i][4]
                    pos_z       = self.mm_config[i][5]              # if chosen atom is 'shel' then get 'core' position

                    rx = self.mm_config[i][3] - self.mm_config[i-1][3]      # r_shel - r_core ... r_cs (r_c_to_s)
                    ry = self.mm_config[i][4] - self.mm_config[i-1][4]
                    rz = self.mm_config[i][5] - self.mm_config[i-1][5]

                    ux = rx*self.mm_config[i][2]
                    uy = ry*self.mm_config[i][2]
                    uz = rz*self.mm_config[i][2]

                    self.atomic_dip_sum[0] += ux
                    self.atomic_dip_sum[1] += uy
                    self.atomic_dip_sum[2] += uz
                    #ux=ux*-1.; uy=uy*-1.; uz=uz*-1.

                    print("%3s%2s%10.6f%12.6f%12.6f%18.6f%14.6f%14.6f" % (atom_name,"s",pos_x,pos_y,pos_z,ux,uy,uz))
                else:
                    atom_name   = self.mm_config[i-1][0]
                    pos_x       = self.mm_config[i-1][3]
                    pos_y       = self.mm_config[i-1][4]
                    pos_z       = self.mm_config[i-1][5]                # if chosen atom is 'shel' then get 'core' position
                    
                    dip_x       = self.mm_dip[i][0] - self.mm_dip[i-1][0]           # Following standard convention ... r(-) - r(+)
                    dip_y       = self.mm_dip[i][1] - self.mm_dip[i-1][1]       # rs*shel - rc*core
                    dip_z       = self.mm_dip[i][2] - self.mm_dip[i-1][2]

                    print("%3s%12.6f%12.6f%12.6f%18.6f%14.6f%14.6f" % (atom_name,pos_x,pos_y,pos_z,dip_x,dip_y,dip_z))
        for i in range(self.qm_cnt):

            if self.sumofatomicdip == 1:

                atom_name       = self.qm_config[i][0]
                pos_x           = self.qm_config[i][3]
                pos_y           = self.qm_config[i][4]
                pos_z           = self.qm_config[i][5]

                ux = 2.*self.mo_config[i][2]*self.mo_config[i][3]*self.pos_integral*self.qm_config[i][2]
                uy = 2.*self.mo_config[i][2]*self.mo_config[i][4]*self.pos_integral*self.qm_config[i][2]
                uz = 2.*self.mo_config[i][2]*self.mo_config[i][5]*self.pos_integral*self.qm_config[i][2]

                self.atomic_dip_sum[0] += ux
                self.atomic_dip_sum[1] += uy
                self.atomic_dip_sum[2] += uz

                #ux=ux*-1; uy=uy*-1.; uz=uz*-1.         

                print("%3s%12.6f%12.6f%12.6f%18.6f%14.6f%14.6f" % (atom_name,pos_x,pos_y,pos_z,ux,uy,uz))

            else:
                atom_name       = self.qm_config[i][0]
                pos_x           = self.qm_config[i][3]
                pos_y           = self.qm_config[i][4]
                pos_z           = self.qm_config[i][5]

                dip_x           = self.qm_dip[i][0] - self.qm_config[i][1]*self.qm_config[i][3]
                dip_y           = self.qm_dip[i][1] - self.qm_config[i][1]*self.qm_config[i][4]
                dip_z           = self.qm_dip[i][2] - self.qm_config[i][1]*self.qm_config[i][5]

                # self.qm_dip[i][x] .EQ. qc*rc + qs*rc + qs*<eigenfunc|rs|eigenfunc>
                # dip_x ... e.g., .EQ. qs*rc + qs*<eigengunc|rs|eigenfunc>

                #print("%3s%12.6f%12.6f%12.6f%18.6e%14.6e%14.6e" % (atom_name,pos_x,pos_y,pos_z,dip_x,dip_y,dip_z))
                print("%3s%12.6f%12.6f%12.6f%18.6f%14.6f%14.6f" % (atom_name,pos_x,pos_y,pos_z,dip_x,dip_y,dip_z))
        print("-"*90)

        # SYSTEM BASIC INFO
        if self.base_if_charge_neutral == False:
            print("")
            print(" @WARNING SYSTEM IS NOT CHARGE NUETRAL, HAVING CHARGE OF %.6e" % (self.base_total_charge_sum))
            print(" CALCULATED DIPOLE MOMEMT MAY VARY BY A CHOICE OF REFERENCE FRAME")
            print("")

        if self.base_if_optimised == False:
            print("")
            print(" @WARNING SYSTEM MAY NOT BE FULLY RELAXED")
            print(" PLEASE CHECK 'Gnorm' or 'Geometric Derivatives' SECTION IN SLAM STANDARD OUTPUT FILE")
            print("")

        # END BASE INFO
        print("")
        print(" Cluster Dipole MM    (e*Angs) :%12.6f%12.6f%12.6f" % (self.cluster_dip_mm[0], self.cluster_dip_mm[1], self.cluster_dip_mm[2]))
        print(" Cluster Dipole QM    (e*Angs) :%12.6f%12.6f%12.6f" % (self.cluster_dip_qm[0], self.cluster_dip_qm[1], self.cluster_dip_qm[2]))
        '''
        if self.sumofatomicdip == 1:
            print(" Total Cluster Dipole (e*Angstrom) :%12.6f%12.6f%12.6f" % (self.atomic_dip_sum[0],self.atomic_dip_sum[1],self.atomic_dip_sum[2]))
        else:
            print(" Total Cluster Dipole (e*Angstrom) :%12.6f%12.6f%12.6f" % (self.cluster_dip[0], self.cluster_dip[1], self.cluster_dip[2]))
        '''
        print(" Total Cluster Dipole (e*Angs) :%12.6f%12.6f%12.6f" % (self.cluster_dip[0], self.cluster_dip[1], self.cluster_dip[2]))
        print("")

        print("#"*90)
        print("%45s" % ("Finalising"))
        print("#"*90)

if __name__=='__main__':

    dip_inst = slam_dipole_mod(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])
    dip_inst.get_base()
    dip_inst.get_cluster_dipole()
    dip_inst.write()
