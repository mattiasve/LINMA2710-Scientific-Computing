
#include "solver.hpp"

/*
    Solve the defined problem using a finite volume method
        with scalar function a, f and g
        with center of the control volume in v_points

        returns the solution sol[i] at points (x[i], y[i])
*/
void solve(const std::function<double(double,double)>& a, const std::function<double(double,double)>& f, const std::function<double(double,double)>& g, std::vector<Point2> v_points, vector<double> &x, vector<double> &y, vector<double> &sol)
{
    int n_points = v_points.size();

    /*
        Compute the finite volumes and the relevant quantities of the neighbours volumes
    */
    Voronoi2 *p_voro; // Voronoi data structure
    std::vector<VoroCell2*> v_voronoi_cells; // list of voronoi cells (=finite volumes)
    Fade_2D *dt; // main pointer to the Delaunay triangulation used by the umericla geometry algorithm, we save it to free the structure later
    compute_volumes(v_points, p_voro, v_voronoi_cells, dt);

    std::vector<std::vector<Neighbour_Volume*>*> adjacent_volumes; // ajacency list of the neighbour volumes, see the Neighbour_Volume class for more information
    treat_volumes( p_voro, v_voronoi_cells, adjacent_volumes);

    /*
        Pose the equations by filling A and b
    */
    arma::SpMat<double> A;
    arma::Col<double> b;
    setup_equations(a, f, g, n_points, v_voronoi_cells, adjacent_volumes, A, b);

    /*
		Solve it
	*/
	arma::Col<double> u = arma::spsolve(A, b, "lapack");


    for (int i=0; i<n_points; i++)
	{
		VoroCell2 *v_cell = v_voronoi_cells[i];
		int idx = v_cell->getCustomCellIndex();
		x.push_back(v_cell->getSite()->x());
		y.push_back(v_cell->getSite()->y());
		sol.push_back(u(idx)); // we use the custom index of the voronoi cell (=finite volume) as the index of the finite volume in the equations
	}

    /*
        Free used data-structures
    */
    delete dt;
    for (size_t i=0; i<adjacent_volumes.size(); i++)
    {
        for (size_t j=0; j<(*adjacent_volumes[i]).size(); j++)
            delete (*adjacent_volumes[i])[j];
        delete adjacent_volumes[i];
    }
}

/*
	Compute the finite volumes with center in v_points using Delaunay triangulation and Voronoi cells
        return the results in p_voro and v_voronoi_cells
*/
void compute_volumes(std::vector<Point2> v_points, Voronoi2* &p_voro, std::vector<VoroCell2*> &v_voronoi_cells, GEOM_FADE2D::Fade_2D* &dt)
{
    dt = new GEOM_FADE2D::Fade_2D();
	dt->insert(v_points);
	p_voro = dt->getVoronoiDiagram();
	p_voro->show("voronoi_and_delaunay.ps",true,false,true,true,false);
	p_voro->show("voronoi_and_sites.ps",true,true,true,false,false); 

	p_voro->getVoronoiCells(v_voronoi_cells);
	for (size_t i=0; i<v_voronoi_cells.size(); i++)
	{
		VoroCell2 *vCell = v_voronoi_cells[i];
		vCell->setCustomCellIndex(i);
	}
}

/*
    Get the neighbours of each volume and compute relevant quantities
        p_voro and v_voronoi_cells are the output of the compute volumes function
        return adjacent volumes an adjacency list of voronoi cell -> (cell index, x, y, distance, common edge length)
*/
void treat_volumes(Voronoi2 *p_voro, std::vector<VoroCell2*> &v_voronoi_cells, std::vector<std::vector<Neighbour_Volume*>*> &adjacent_volumes)
{
    for (size_t i=0; i<v_voronoi_cells.size(); i++)
    {
        VoroCell2 *v_cell = v_voronoi_cells[i];
        // std::cout<<v_cell->getCustomCellIndex()<<std::endl;
        // Point2 *site = v_cell->getSite();
        // std::cout<<*site<<std::endl;
        // int area = v_cell->getArea(); -----> AREA 
        // std::cout<<"\tcell area = ";
        // if (area>0.)
        //     std::cout<<area<<std::endl;
        // else
        //     std::cout<<"+Inf (unbounded cell)"<<std::endl;


        // std::cout<<"\tadjacent voronoi vertices : "<<std::endl;
        std::vector<VoroVertex2*> v_voro_vertices;
        v_cell->getVoronoiVertices(v_voro_vertices);
        // for (size_t j=0; j<v_voro_vertices.size(); j++)
        // {
        //     VoroVertex2 *v_vertex = v_voro_vertices[j];
        //     std::cout<<"\t\t"<<v_vertex->getPoint()<<std::endl;
        // }

        // std::cout<<"\tvoronoi neighbour sites with non-zero edge length :"<<std::endl;
        adjacent_volumes.push_back(new std::vector<Neighbour_Volume*>());
        for (size_t j=0; j<v_voro_vertices.size(); j++)
        {
            VoroCell2* p_cell0;
            VoroCell2* p_cell1;
            Point2 p1 = v_voro_vertices[j]->getPoint();
            Point2 p2 = v_voro_vertices[(j+1)%v_voro_vertices.size()]->getPoint();
            double length = distance(&p1, &p2);
            if (length>1e-10 && p_voro->getVCellsAtVEdge(v_voro_vertices[j], v_voro_vertices[(j+1)%v_voro_vertices.size()], p_cell0, p_cell1))
            {
                double dist0 = distance(p_cell0->getSite(), v_cell->getSite());
                double dist1 = distance(p_cell1->getSite(), v_cell->getSite());
                double dist;
                VoroCell2* p_cell;
                if (dist0 > 1e-10)
                {
                    p_cell = p_cell0;
                    dist = dist0;
                }
                if (dist1 > 1e-10)
                {
                    p_cell = p_cell1;
                    dist = dist1;
                }
                Point2 *point = p_cell->getSite();
                adjacent_volumes[i]->push_back(new Neighbour_Volume(p_cell->getCustomCellIndex(), (v_cell->getSite()->x()+point->x())/2, (v_cell->getSite()->y()+point->y())/2 ,dist, length));
                // std::cout<<"\t\t"<<p_cell->getCustomCellIndex()<<"\t"<<*(p_cell->getSite())<<std::endl;
                // std::cout<<"\t\t\tdistance with neighbour h_i,j = "<<dist<<std::endl;
                // std::cout<<"\t\t\tcommon length with neighbour l_i,j = "<<length<<std::endl;
            }
        }
        // std::cout<<std::endl;
    }
}

/*
    Setup system of equations
*/
void setup_equations(const std::function<double(double,double)>& a, const std::function<double(double,double)>& f, const std::function<double(double,double)>& g, int n_points, std::vector<VoroCell2*> &v_voronoi_cells, std::vector<std::vector<Neighbour_Volume*>*> &adjacent_volumes, arma::SpMat<double> &A, arma::Col<double> &b)
{
	// TODO
    A = arma::sp_mat(n_points, n_points);
    b = arma::vec(n_points);
    // A.print();
    std::cout << "length of voronoi cell list = " << v_voronoi_cells.size() << std::endl;

    for (int i=0; i<n_points; i++)
    {
        int nb_neigbour = (*adjacent_volumes[i]).size();
        
        if (v_voronoi_cells[i]->getArea() == -1) // cell is on the boundary, apply boundary condition
        {
            A(i, v_voronoi_cells[i]->getCustomCellIndex()) = 1; 
            b[i] = g(v_voronoi_cells[i]->getSite()->x(), v_voronoi_cells[i]->getSite()->y());
        }
        else {
            for (int j=0; j<nb_neigbour; j++) // cell is inside the domain
            {
                double Aij;
                double xij = (*adjacent_volumes[i])[j]->x;
                double yij = (*adjacent_volumes[i])[j]->y;
                // compute coefficients and fill matrix
                Aij = a(xij, yij) * ((*adjacent_volumes[i])[j]->l)/(*adjacent_volumes[i])[j]->h ;
                A(i, (*adjacent_volumes[i])[j]->idx) += Aij;
                A(i, v_voronoi_cells[i]->getCustomCellIndex()) -= Aij;
                // compute right-hand side and fill dense vector
                b[i] = f(xij, yij) * v_voronoi_cells[i]->getArea();
            }
        }
    }
    //std::cout << "IM HERE" << std::endl;
    //A.print();
}

/*
    computes the classical euclidean distance between two points
*/
double distance(Point2 *p0, Point2 *p1)
{
	return sqrt(pow((p0->x()-p1->x()),2) + pow((p0->y()-p1->y()),2));
}