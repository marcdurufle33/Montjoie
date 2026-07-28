from numpy import *
from quadrature import *

class HelmholtzSolver:
    """ Classe pour calculer les matrices elements finis pour 
    un rectangle (maillage regulier) """
    def __init__(self, Nx, Ny, xmin, xmax, ymin, ymax, r):
        self.pos_x = linspace(xmin, xmax, Nx+1)
        self.pos_y = linspace(ymin, ymax, Ny+1)
        self.Nodle = self.CreateRegularNumbering(Nx, Ny, r)
        self.nodl = (Nx*r+1)*(Ny*r+1)
        self.points_lob, self.poids_lob = ComputeLobattoJacobi(r, 0, 0)
        self.points_gauss, self.poids_gauss = ComputeGaussJacobi(r-1, 0, 0)
        print("points Gauss-Lobatto = ", self.points_lob)
        self.phi = LagrangeFunctions(self.points_lob)
        #for i in range (Nx):
        #    for j in range(Ny):
        #        print ("Numerotation", i, j, self.Nodle[j*Nx+i])
    
    def CreateRegularNumbering(self, Nx, Ny, r):
        """ Creation de la numerotation pour chaque element
        cas d'un maillage regulier d'un rectangle """
        Nodle = [0]*Nx*Ny
        for i in range (Nx):
            for j in range(Ny):
                num = j*Nx + i
                arr = zeros((r+1)*(r+1), dtype='int64')
                # 4 sommets
                arr[0] = j*(Nx+1)+i
                arr[r] = j*(Nx+1)+i+1
                arr[r*(r+1)] = (j+1)*(Nx+1)+i
                arr[r*(r+1)+r] = (j+1)*(Nx+1)+i+1
                
                # ddls des aretes
                off_edge_x = (Nx+1)*(Ny+1) + (j*Nx + i)*(r-1)
                off_edge_y = (Nx+1)*(Ny+1) + Nx*(Ny+1)*(r-1) 
                off_edge_y += (j*(Nx+1)+i)*(r-1)
                for k in range(1, r):
                    arr[k] = off_edge_x + k-1
                    arr[k*(r+1)] = off_edge_y + k-1
                    arr[k*(r+1)+r] = off_edge_y + r-1 +k-1
                    arr[r*(r+1) + k] = off_edge_x + Nx*(r-1) + k-1
                    
                # ddls a l'interieur
                off_elt = (Nx+1)*(Ny+1) + Nx*(Ny+1)*(r-1)+ Ny*(Nx+1)*(r-1) 
                off_elt += num*(r-1)*(r-1)
                for k in range(1, r):
                    for m in range(1, r):
                        arr[k*(r+1) + m] = off_elt + (k-1)*(r-1) + m-1
                        
                Nodle[num] = arr
                
        return Nodle
    
    def ComputeStiffnessMatrix(self):
        """ Calcul de la matrice rigidite (comme une matrice dense) """
        K = zeros([self.nodl, self.nodl])
        r = len(self.points_lob)-1
        Nx = len(self.pos_x)-1
        Ny = len(self.pos_y)-1
        grad = self.phi.ComputeGradPhi()
        #print ("gradient = ", grad)
        Kelem = zeros([r+1, r+1])
        for i in range(r+1):
            for j in range(r+1):
                for k in range(r+1):
                    Kelem[i, j] += self.poids_lob[k]*grad[i, k]*grad[j, k]
        
        print ("matrice rigidite 1-D = ", Kelem)
        C0 = 1.0; C1 = 1.0; # cas isotrope
        for ie in range(Nx):
            for je in range(Ny):
                num = je*Nx + ie
                dx = self.pos_x[ie+1] - self.pos_x[ie]
                dy = self.pos_y[je+1] - self.pos_y[je]
                C0h = C0*dy / dx; C1h = C1 * dx / dy;
                for i1 in range(r+1):
                    for i2 in range(r+1):
                        i = i2*(r+1) + i1
                        for j1 in range(r+1):
                            j = i2*(r+1)+j1
                            #print("numeros", num, i, j, len(self.Nodle))
                            K[self.Nodle[num][i], self.Nodle[num][j]] += C0h*Kelem[i1, j1]*self.poids_lob[i2]
                            
                        for j2 in range(r+1):
                            j = j2*(r+1) + i1;
                            K[self.Nodle[num][i], self.Nodle[num][j]] += C1h*Kelem[i2, j2]*self.poids_lob[i1]
        return K

    def ComputeMassMatrix(self):
        """ Calcul de la matrice de masse (diagonale) """
        M = zeros(self.nodl)
        r = len(self.points_lob)-1
        Nx = len(self.pos_x)-1
        Ny = len(self.pos_y)-1
        for ie in range(Nx):
            for je in range(Ny):
                dx = self.pos_x[ie+1] - self.pos_x[ie]
                dy = self.pos_y[je+1] - self.pos_y[je]
                num = je*Nx + ie
                jacob = dx*dy
                for i1 in range(r+1):
                    for i2 in range(r+1):                        
                        i = i2*(r+1) + i1
                        M[self.Nodle[num][i]] += jacob*self.poids_lob[i1]*self.poids_lob[i2]
                        
        return M
        
    def InterpolateField(self, XI, YI, V):
        """ Interpolation d'un champ sur des points 2-D """
        Vi = zeros(XI.shape)
        # cas regulier
        dx = self.pos_x[1] - self.pos_x[0]
        dy = self.pos_y[1] - self.pos_y[0]
        Nx = len(self.pos_x)-1
        Ny = len(self.pos_y)-1
        r = len(self.points_lob)-1
        for i in range(XI.shape[0]):
            for j in range(XI.shape[1]):
                x = XI[i, j]; y = YI[i, j]
                ie = int(floor((x - self.pos_x[0])/dx))
                je = int(floor((y - self.pos_y[0])/dy))
                if ((ie > 0) and (ie <= Nx) and (x == self.pos_x[ie])):
                    ie -= 1
                    
                if ((je > 0) and (je <= Ny) and (y == self.pos_y[je])):
                    je -= 1
                
                if ((ie >= 0) and (ie < Nx) and (je >= 0) and (je < Ny)):
                    num = je*Nx + ie
                    xloc = (x - self.pos_x[ie])/dx
                    yloc = (y - self.pos_y[je])/dy
                    #print("x, y", x, y, ie, je, xloc, yloc)
                    phi_x = self.phi.ComputeValuesPhiRef(xloc)
                    phi_y = self.phi.ComputeValuesPhiRef(yloc)
                    for k1 in range(r+1):
                        for k2 in range(r+1):
                            k = k2*(r+1) + k1
                            ndof = self.Nodle[num][k]
                            Vi[i, j] += phi_x[k1]*phi_y[k2]*V[ndof]
        
        return Vi
        
    def ComputeFirstOrderSystem(self):
        """ calcul du systeme du premier ordre
        en prenant les points de Gauss pour l' inconnue vectorielle, Gauss-Lobatto pour l'inconnue scalaire """
        Nx = len(self.pos_x)-1
        Ny = len(self.pos_y)-1
        r = len(self.points_lob)-1
        nodl_vec = Nx*Ny*r*r        
        nb_dof = self.nodl + 2*nodl_vec
        M = zeros(nb_dof)
        K = zeros([nb_dof, nb_dof])
        
        # gradient des fcts de base sur les points de Gauss
        grad = zeros([r+1, r])
        val_phiG = zeros([r+1, r])
        grad_GL = self.phi.ComputeGradPhi()
        for k in range(r):
            eval_phi = self.phi.ComputeValuesPhiRef(self.points_gauss[k])
            for i in range(r+1):
                val_phiG[i, k] = eval_phi[i]
                for j in range(r+1):
                    grad[i, k] += grad_GL[i, j] * eval_phi[j]
        
        # boucle sur les elements
        for ie in range(Nx):
            for je in range(Ny):
                num = je*Nx + ie
                dx = self.pos_x[ie+1] - self.pos_x[ie]
                dy = self.pos_y[je+1] - self.pos_y[je]
                jacob = dx*dy
                for i1 in range(r+1):
                    for i2 in range(r+1):
                        i = i2*(r+1) + i1
                        ig = self.Nodle[num][i]
                        # matrice de masse scalaire
                        M[ig] += jacob*self.poids_lob[i1]*self.poids_lob[i2]
                        # matrice de rigidite
                        for j1 in range(r):
                            for j2 in range(r):
                                j = j2*r + j1
                                jg = j + num*r*r
                                poids = self.poids_gauss[j1]*self.poids_gauss[j2]
                                K[ig, self.nodl + 2*jg] = dy*grad[i1, j1]*val_phiG[i2, j2]*poids
                                K[ig, self.nodl + 2*jg+1] = dx*grad[i2, j2]*val_phiG[i1, j1]*poids
                                K[self.nodl + 2*jg, ig] = -K[ig, self.nodl + 2*jg]
                                K[self.nodl + 2*jg+1, ig] = -K[ig, self.nodl + 2*jg+1]
                
                for i1 in range(r):
                    for i2 in range(r):
                        i = i2*r + i1
                        ig = i + num*r*r
                        # matrice de masse vectorielle
                        M[self.nodl + 2*ig] += jacob*self.poids_gauss[i1]*self.poids_gauss[i2]
                        M[self.nodl + 2*ig+1] += jacob*self.poids_gauss[i1]*self.poids_gauss[i2]
        
        return K, M

