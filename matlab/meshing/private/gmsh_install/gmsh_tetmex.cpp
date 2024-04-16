#include "mex.hpp"
#include "mexAdapter.hpp"
#include <stdio.h>
#include <iostream>
#include <vector>

#include "gmsh.h"

using namespace matlab::data;
using namespace std;

#define in_bb(bounds,x,y,z)(!(x<bounds[0] || x>bounds[3] || y<bounds[1] || y>bounds[4] || z<bounds[2] || z>bounds[5]));

// class for the mexfunction so that it may be called
class MexFunction : public matlab::mex::Function {
    // Pointer to MATLAB engine to call fprintf
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    // Factory to create MATLAB data arrays
    ArrayFactory factory;

    // Create an output stream
    std::ostringstream stream;
public:
     // operator so function may be called in matlab as tetgen_mex()
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        // getting number of arguments
        size_t nargin = inputs.size();
        size_t nargout = outputs.size();

        // grabbing model string
        std::string modelname = matlab::engine::convertUTF16StringToUTF8String(inputs[1][0]);
        
        // grabbing average edge length
        int n = std::move(inputs[0][0]);
        double h = 1/((double)n-1);

        // grabbing min edge ratio for h-GR
        double h_ratio = std::move(inputs[2][0]);
        

        // calling gmsh to create surface mesh
        call_gmsh_volume(modelname, h, h_ratio, outputs);

        
    }

    void call_gmsh_volume(const std::string& modelname, const double &h, const double &h_ratio, matlab::mex::ArgumentList &outputs){

        // finding h-gr targets
        vector<double> target_coords, h_targets;
        find_meshsize_based_curvature_3D(target_coords, h_targets, modelname, h_ratio*h, h, 15.0);

        gmsh::initialize();
        gmsh::model::add("mesh_"+modelname);
        gmsh::merge(modelname+".geo");
        gmsh::vectorpair entities;
        gmsh::model::getEntities(entities, 3);

        // adding box points and removing volume
        int vs[8];
        vs[0] = gmsh::model::geo::addPoint(0.0,0.0,0.0);
        vs[1] = gmsh::model::geo::addPoint(1.0,0.0,0.0);
        vs[2] = gmsh::model::geo::addPoint(1.0,1.0,0.0);
        vs[3] = gmsh::model::geo::addPoint(0.0,1.0,0.0);
        vs[4] = gmsh::model::geo::addPoint(0.0,0.0,1.0);
        vs[5] = gmsh::model::geo::addPoint(1.0,0.0,1.0);
        vs[6] = gmsh::model::geo::addPoint(1.0,1.0,1.0);
        vs[7] = gmsh::model::geo::addPoint(0.0,1.0,1.0);
        gmsh::model::geo::synchronize();

        int E[12];
        E[0] = gmsh::model::geo::addLine(vs[0],vs[1]);
        E[1] = gmsh::model::geo::addLine(vs[1],vs[2]);
        E[2] = gmsh::model::geo::addLine(vs[2],vs[3]);
        E[3] = gmsh::model::geo::addLine(vs[3],vs[0]);
        E[4] = gmsh::model::geo::addLine(vs[0],vs[4]);
        E[5] = gmsh::model::geo::addLine(vs[1],vs[5]);
        E[6] = gmsh::model::geo::addLine(vs[2],vs[6]);
        E[7] = gmsh::model::geo::addLine(vs[3],vs[7]);
        E[8] = gmsh::model::geo::addLine(vs[4],vs[5]);
        E[9] = gmsh::model::geo::addLine(vs[5],vs[6]);
        E[10] = gmsh::model::geo::addLine(vs[6],vs[7]);
        E[11] = gmsh::model::geo::addLine(vs[7],vs[4]);
        gmsh::model::geo::synchronize();

        int S[6],s;
        s = gmsh::model::geo::addCurveLoop({-E[3], -E[2], -E[1], -E[0]},-1,true);
        S[0] = gmsh::model::geo::addPlaneSurface({s});
        s = gmsh::model::geo::addCurveLoop({E[0], E[5], -E[8], -E[4]},-1,true);
        S[1] = gmsh::model::geo::addPlaneSurface({s});
        s = gmsh::model::geo::addCurveLoop({E[1], E[6], -E[9], -E[5]},-1,true);
        S[2] = gmsh::model::geo::addPlaneSurface({s});
        s = gmsh::model::geo::addCurveLoop({E[2], E[7], -E[10], -E[6]},-1,true);
        S[3] = gmsh::model::geo::addPlaneSurface({s});
        s = gmsh::model::geo::addCurveLoop({E[3], E[4], -E[11], -E[7]},-1,true);
        S[4] = gmsh::model::geo::addPlaneSurface({s});
        s = gmsh::model::geo::addCurveLoop({E[8], E[9], E[10], E[11]},-1,true);
        S[5] = gmsh::model::geo::addPlaneSurface({s});
        gmsh::model::geo::synchronize();

        s = gmsh::model::geo::addSurfaceLoop({S[0],S[1],S[2],S[3],S[4],S[5]});
        vector<int> sltags;
        sltags.push_back(s);
        for (auto e : entities){
            sltags.push_back(e.second);
        }
        gmsh::model::geo::addVolume(sltags);
        gmsh::model::geo::synchronize();
        gmsh::model::geo::remove(entities);
        gmsh::model::geo::synchronize();

        // setting value
        gmsh::option::setNumber("General.Verbosity", 3);
        gmsh::option::setNumber("Mesh.MeshSizeMax",h);
        gmsh::option::setNumber("Mesh.MeshSizeMin",h_ratio*h);
        gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary",0);
        gmsh::option::setNumber("Mesh.MeshSizeFromPoints",0);
        gmsh::option::setNumber("Mesh.MeshSizeFromCurvature",0);

        create_curvature_bgmesh(target_coords, h_targets, h, h_ratio*h);

        gmsh::model::mesh::generate(3);
        gmsh::model::mesh::removeDuplicateNodes();
        stream << "finished generating 3d mesh" << endl;
        displayOnMATLAB(stream);

        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t> > elementTags;
        std::vector<std::vector<std::size_t> > nodeTags_Elems;
        std::vector<std::size_t>  nodeTags;
        std::vector<double> coord;
        std::vector<double> parametricCoord;
        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags_Elems, 3);
        gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, 3,-1, true,false);
        
        int nv = nodeTags.size();
        
        vector<double> coords(3*nv);
        int np = 0;
        int maxn = 0;
        for (int n = 0; n<nodeTags.size(); ++n){
            if (nodeTags[n] > maxn){
                maxn = nodeTags[n];
            }
        }
        maxn = maxn+10;
        vector<int> map(maxn,-1);
        for (int n = 0; n<nodeTags.size(); ++n){
            if (map[nodeTags[n]]==-1){
                map[nodeTags[n]] = np;
                for (int j = 0;j<3;++j){
                    coords[3*np+j] = coord[3*n+j];
                }
                ++np;
            }
        }

        TypedArray<double> xs = factory.createArray<double>({(long long unsigned int)np, 3});
        for (int n = 0; n<np; ++n){
            for (int j = 0; j<3; ++j){
                xs[n][j] = coords[3*n+j];
            }
        }

        std::string elemname;
        int dim, order, numNodes, numPrimaryNodes, nelems, type;
        std::vector<double> localNodeCoord;
        gmsh::model::mesh::getElementProperties(elementTypes[0], elemname, dim, order, numNodes, localNodeCoord, numPrimaryNodes);
        nelems = nodeTags_Elems[0].size()/numNodes;

        TypedArray<int32_t> elems = factory.createArray<int32_t>({(long long unsigned int)nelems, (long long unsigned int)numNodes});
        for (int e = 0; e<nelems; ++e){
            for (int n = 0; n<numNodes; ++n){
                    elems[e][n] = map[nodeTags_Elems[0][numNodes*e + n]]+1;
            }
        }


        gmsh::finalize();

        outputs[0] = elems;
        outputs[1] = xs;
    }

    void find_meshsize_based_curvature_3D(vector<double> &target_coords, vector<double> &h_targets, string modelname, double min_sz, double max_sz, double theta){
        // creating a mesh of the hole
        gmsh::initialize();
        gmsh::option::setNumber("General.Verbosity", 3);
        gmsh::model::add("mesh_"+modelname);
        gmsh::merge(modelname+".geo");
        gmsh::option::setNumber("Mesh.MeshSizeMax",max_sz);
        gmsh::option::setNumber("Mesh.MeshSizeMin",min_sz);
        gmsh::model::mesh::generate(2);

        // resizing arrays
        target_coords.clear();
        h_targets.clear();

        // getting principal curvatures
        double h_target;
        gmsh::vectorpair entities;
        gmsh::model::getEntities(entities, 2);
        vector<size_t> tags;
        vector<double> coord,param,CurvMax, CurvMin, DirMax, DirMin;
        for (auto e : entities){
            int surftag = e.second;

            gmsh::model::mesh::getNodes(tags, coord, param, 2, surftag, true);
            gmsh::model::getPrincipalCurvatures(surftag, param, CurvMax, CurvMin, DirMax, DirMin);

            if (modelname.find("sphere") != std::string::npos){
                for (int i = 0; i<CurvMax.size(); ++i){
                    CurvMax[i] = 10.0;
                }
            }

            int nv = coord.size()/3;
            for (int i = 0; i<nv; ++i){
                h_target = min(max_sz, max(min_sz, (theta*M_PI/180.0) / abs(CurvMax[i])));
                if (h_target != max_sz){
                    target_coords.push_back(coord[3*i]);
                    target_coords.push_back(coord[3*i+1]);
                    target_coords.push_back(coord[3*i+2]);
                    h_targets.push_back(h_target);
                }
            }
        }

        gmsh::finalize();
    }

    void create_curvature_bgmesh(const vector<double> &target_coords, const vector<double> &h_targets, double max_sz, double min_sz){
        int ngroups = (int)round(max_sz/min_sz);
        double dg = max_sz/min_sz;
        // no bgmesh if no refinement needed
        if (ngroups == 1){return;}

        vector<vector<double>> groups(ngroups-1);

        int npoints = h_targets.size(), tag_node, tag_dist, tag_thres;
        vector<double> Fieldlist;
        int np = 0,group;
        for (int i = 0; i<npoints; ++i){
            group = (int)(round(h_targets[i]/min_sz)-1);
            if (group < ngroups-1){
                tag_node = gmsh::model::geo::addPoint(target_coords[3*i],target_coords[3*i+1],target_coords[3*i+2]);
                groups[group].push_back(tag_node);
            }
        }

        double n = 1;
        for (auto g : groups){
            // add distance field
            if (g.size() > 0){
            tag_dist = gmsh::model::mesh::field::add("Distance");
            gmsh::model::mesh::field::setNumbers(tag_dist, "PointsList", g);

            // add threshold field
            tag_thres = gmsh::model::mesh::field::add("Threshold");
            gmsh::model::mesh::field::setNumber(tag_thres, "InField", tag_dist);
            gmsh::model::mesh::field::setNumber(tag_thres, "SizeMin", n*min_sz);
            gmsh::model::mesh::field::setNumber(tag_thres, "SizeMax", max_sz);
            gmsh::model::mesh::field::setNumber(tag_thres, "DistMin", 0.0);
            gmsh::model::mesh::field::setNumber(tag_thres, "DistMax", max_sz);

            Fieldlist.push_back(tag_thres);
        
            ++n;
            }
        }

        int bgtag = gmsh::model::mesh::field::add("Min");
        gmsh::model::mesh::field::setNumbers(bgtag, "FieldsList", Fieldlist);

        int octree = gmsh::model::mesh::field::add("Octree");
        gmsh::model::mesh::field::setNumber(octree, "InField", (double)bgtag);
        gmsh::model::mesh::field::setAsBackgroundMesh(octree);
        gmsh::model::geo::synchronize();
    }



    void displayOnMATLAB(std::ostringstream& stream) {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
            std::vector<Array>({ factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }
};