/*****************************************************************************
*
* Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/
/* SIMPLE SIMULATION SKELETON */
#include <VisItControlInterface_V2.h>
#include <VisItDataInterface_V2.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//#define DEBUG_PRINT

#include <iostream>
#include <sstream>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "patch.h"
#include "SimulationExample.h"

#define VISIT_COMMAND_PROCESS 0
#define VISIT_COMMAND_SUCCESS 1
#define VISIT_COMMAND_FAILURE 2

/* Data Access Function prototypes */
visit_handle SimGetMetaData(void *);
visit_handle SimGetMesh(int, const char *, void *);
visit_handle SimGetVariable(int, const char *, void *);
visit_handle SimGetDomainNesting(const char *, void *);
visit_handle SimGetDomainList(const char *, void *);

/******************************************************************************
 * Code for calculating data values
 ******************************************************************************/

class complex
{
public:
    complex() : a(0.), b(0.) { }
    complex(float A, float B) : a(A), b(B) { }
    complex(const complex &obj) : a(obj.a), b(obj.b) { }
    complex operator = (const complex &obj) { a = obj.a; b = obj.b; return *this;}
    complex operator + (const complex &obj) const
    {
        return complex(a + obj.a,  b + obj.b);
    }
    complex operator * (const complex &obj) const
    {
        return complex(a * obj.a - b * obj.b, a * obj.b + b * obj.a);
    }
    float mag2() const
    {
        return a*a + b*b;
    }
    float mag() const
    {
        return sqrt(a*a + b*b);
    }
private:
    float a,b;
};

#define MAXIT 30

int
mandelbrot(const complex &C)
{
    complex Z;
    for(int zit = 0; zit < MAXIT; ++zit)
    {
        Z = (Z * Z) + C;
        if(Z.mag2() > 4.f)
            return zit+1;
    }
    return 0;
}

void
patch_calculate_data(patch_t *patch)
{
    int i,j;
    float ty,y,tx,x;
    float cellWidth, x0, x1;
    float cellHeight, y0, y1;
    unsigned char *data = patch->data;
    
    /* Compute x0, x1 and y0,y1 which help us locate cell centers. */
    cellWidth = (patch->window[1] - patch->window[0]) / ((float)patch->nx);
    x0 = patch->window[0] + cellWidth / 2.f;
    x1 = patch->window[1] - cellWidth / 2.f;
    cellHeight = (patch->window[3] - patch->window[2]) / ((float)patch->ny);
    y0 = patch->window[2] + cellHeight / 2.f;
    y1 = patch->window[3] - cellHeight / 2.f;
    for(j = 0; j < patch->ny; ++j)
    {
        ty = (float)j / (float)(patch->ny - 1);
        y = (1.f-ty)*y0 + ty*y1;
        for(i = 0; i < patch->nx; ++i)
        {
            tx = (float)i / (float)(patch->nx - 1);
            x = (1.f-tx)*x0 + tx*x1;

            *data++ = (unsigned char)mandelbrot(complex(x, y));
        }
    }
}

/******************************************************************************
 * Code for helping calculate AMR refinement
 ******************************************************************************/

int
neighbors(patch_t *patch, int i, int j)
{
    const float kernel[3][3] = {
    {0.08f, 0.17f, 0.08f},
    {0.17f, 0.f,   0.17f},
    {0.08f, 0.17f, 0.08f}
    };

    float sum = 0;
    for(int jj = 0; jj < 3; ++jj)
    {
        int J = j + jj - 1;
        for(int ii = 0; ii < 3; ++ii)
        {
            int I = i + ii - 1;
            float value = (float)patch->data[J*patch->nx+I];
            sum += value * kernel[jj][ii];
        }
    }
    return (int)sum;
}

void
detect_refinement(patch_t *patch, image_t *mask)
{
    int i,j;
    /* Let's look for large differences within a kernel. This lets us
     * figure out areas that we need to refine because they contain
     * features. We set a 1 into the mask for cells that need refinement
     */
    for(j = 1; j < patch->ny-1; ++j)
        for(i = 1; i < patch->nx-1; ++i)
        {
            int index = j*patch->nx+i;
            int dval = (int)patch->data[index] - neighbors(patch, i, j);
            if(dval < 0) dval = -dval;
            if(dval > 2) 
                mask->data[index] = 1;
            else
                mask->data[index] = 0;
        }
}

/******************************************************************************
 * Calculate the Mandelbrot set using AMR
 ******************************************************************************/

void
calculate_amr_helper(patch_t *patch, int level, int max_levels, int ratio, int *patchid)
{
    int i;

    patch->id = *patchid;
    *patchid = *patchid + 1;
    patch->level = level;

    /* Calculate the data on this patch */
    patch_calculate_data(patch);

    if(level+1 > max_levels)
        return;

    /* Examine this patch's data and see if we need to refine. If we
     * do refine then we'll create subpatches that we can operate on.
     */
    patch_refine(patch, ratio, detect_refinement);

    for(i = 0; i < patch->nsubpatches; ++i)
        calculate_amr_helper(&patch->subpatches[i], level+1, max_levels, ratio, patchid);
}

void
calculate_amr(patch_t *patch, int level, int max_levels, int ratio)
{
    int patchid = 0;
    calculate_amr_helper(patch, level, max_levels, ratio, &patchid);
}

/******************************************************************************
 * Simulation data and functions
 ******************************************************************************/

#define SIM_STOPPED       0
#define SIM_RUNNING       1

enum {ASSIGNMENT_STATIC, ASSIGNMENT_DYNAMIC, ASSIGNMENT_RANDOM, ASSIGNMENT_RANK0};

struct simulation_data
{
    int     par_rank;
    int     par_size;
    int     cycle;
    double  time;
    int     runMode;
    int     done;
    int     max_levels;
    int     refinement_ratio;
    int     savingFiles;
    int     saveCounter;
    bool    autoupdate;
    bool    echo;
    patch_t patch;

    int     *patch_list; /* Patch list for this rank. */
    int      npatch_list;
    int      patch_assignment;
    int      patch_verbose;
};

void
simulation_data_ctor(simulation_data *sim)
{
    sim->cycle = 0;
    sim->time = 0.;
    sim->runMode = SIM_STOPPED;
    sim->done = 0;
    sim->max_levels = 2;
    sim->refinement_ratio = 4;
    sim->savingFiles = 0;
    sim->saveCounter = 0;
    sim->autoupdate = false;
    sim->echo = false;
    patch_ctor(&sim->patch);

    sim->patch_list = NULL;
    sim->npatch_list = 0;
    sim->patch_assignment = ASSIGNMENT_STATIC;
    sim->patch_verbose = 0;
}

void
simulation_data_dtor(simulation_data *sim)
{
    patch_dtor(&sim->patch);
    if(sim->patch_list != NULL)
    {
        free(sim->patch_list);
        sim->patch_list = NULL;
    }
    sim->npatch_list = 0;
}

const char *cmd_names[] = {"halt", "step", "run", "update", "toggleupdates"};

#ifdef PARALLEL
/******************************************************************************
 *
 * Purpose: Fill in the patch_list so we know which ranks own which patches.
 *
 * Programmer: Brad Whitlock
 * Date:       Mon Aug  6 15:15:35 PDT 2012
 *
 * Modifications:
 *
 *****************************************************************************/

void
map_patches_to_processors(simulation_data *sim)
{
    /* Divide the patches among processors. */
    patch_t **patches = patch_flat_array(&sim->patch);
    int npatches = patch_num_patches(&sim->patch);

    FREE(sim->patch_list);
    sim->patch_list = (int *)malloc(sizeof(int) * npatches);
    memset(sim->patch_list, 0, sizeof(int) * npatches);
    sim->npatch_list = 0;
    if(sim->patch_assignment == ASSIGNMENT_RANK0)
    {
        for(int i = 0; i < npatches; ++i)
        {
            /* All on rank 0. */
            if(sim->par_rank == 0)
            {
                sim->patch_list[sim->npatch_list++] = patches[i]->id;
            }
        }
    }
    else if(sim->patch_assignment == ASSIGNMENT_STATIC)
    {
        for(int i = 0; i < npatches; ++i)
        {
            /* Static assignment based on rank. */
            if(patches[i]->id % sim->par_size == sim->par_rank)
            {
                sim->patch_list[sim->npatch_list++] = patches[i]->id;
            }
        }
    }
    else if(sim->patch_assignment == ASSIGNMENT_DYNAMIC)
    {
        for(int i = 0; i < npatches; ++i)
        {
            /* Pathologically shuffle domains to different processors based 
             * on the cycle to make sure we can handle that.
             */
            if((sim->cycle + patches[i]->id) % sim->par_size == sim->par_rank)
            {
                sim->patch_list[sim->npatch_list++] = patches[i]->id;
            }
        }
    }
    else if(sim->patch_assignment == ASSIGNMENT_RANDOM)
    {
        /* Get the random work assignment from rank 0. */
        int *iptr = (int *)malloc(sizeof(int) * npatches);
        for(int i = 0; i < npatches; ++i)
            iptr[i] = rand() % sim->par_size;
        MPI_Bcast(iptr, npatches, MPI_INT, 0, MPI_COMM_WORLD);

        sim->npatch_list = 0;
        FREE(sim->patch_list);
        sim->patch_list = (int *)malloc(sizeof(int) * npatches);
        for(int i = 0; i < npatches; ++i)
        {
            if(iptr[i] == sim->par_rank)
                sim->patch_list[sim->npatch_list++] = i;
        }
        FREE(iptr);
    }
    FREE(patches);

    if(sim->patch_verbose)
    {
        std::ostringstream s;
        s << "Rank " << sim->par_rank << " domain list: ";
        for(int i = 0; i < sim->npatch_list; ++i)
            s << ", " << sim->patch_list[i];
        std::cout << s.str() << std::endl;
    }
}
#endif

/******************************************************************************
 *
 * Purpose: Called to simulate one timestep of data.
 *
 * Programmer: Brad Whitlock
 * Date:       Thu Mar 19 12:14:23 PDT 2009
 *
 * Modifications:
 *
 *****************************************************************************/

void
simulate_one_timestep(simulation_data *sim)
{
    if(sim->par_rank == 0)
        printf("Simulating time step: cycle=%d, time=%lg\n", sim->cycle, sim->time);

    const float window0[] = {-1.6f, 0.6f, -1.1f, 1.1f};
#define ORIGINX -1.5f
#define ORIGINY -0.5f
#define WSIZE 0.5f
    const float window1[] = {ORIGINX, ORIGINX + WSIZE, ORIGINY, ORIGINY + WSIZE};
#define NX 256
#define NY 256

    /* oscillate between 2 windows */
    float window[4];
    float t = 0.5f * sin(sim->time) + 0.5f;
    window[0] = (1.f - t)*window0[0] + t*window1[0];
    window[1] = (1.f - t)*window0[1] + t*window1[1];
    window[2] = (1.f - t)*window0[2] + t*window1[2];
    window[3] = (1.f - t)*window0[3] + t*window1[3];

    /* Blow away the previous patch data and calculate. */
    patch_dtor(&sim->patch);
    patch_ctor(&sim->patch);
    sim->patch.window[0] = window[0];
    sim->patch.window[1] = window[1];
    sim->patch.window[2] = window[2];
    sim->patch.window[3] = window[3];
    sim->patch.logical_extents[0] = 0;
    sim->patch.logical_extents[1] = NX-1;
    sim->patch.logical_extents[2] = 0;
    sim->patch.logical_extents[3] = NY-1;
    patch_alloc_data(&sim->patch, NX, NY);
    calculate_amr(&sim->patch, 0, sim->max_levels, sim->refinement_ratio);

#ifdef PARALLEL
    map_patches_to_processors(sim);
#endif

    ++sim->cycle;
    sim->time += (M_PI / 30.);
}


void read_input_deck(void) { }


/******************************************************************************
 *
 * Purpose: This is the main event loop function.
 *
 * Programmer: Brad Whitlock
 * Date:       Thu Mar 19 12:14:23 PDT 2009
 *
 * Modifications:
 *  Cyrus Harrison,Fri Feb  2 14:06:33 PST 2018
 *  Modify to change into batch example
 * 
 *****************************************************************************/

void mainloop_batch(simulation_data *sim)
{
    /* Explicitly load VisIt runtime functions and install callbacks. */
    VisItInitializeRuntime();
    VisItSetGetMetaData(SimGetMetaData, (void*)sim);
    VisItSetGetMesh(SimGetMesh, (void*)sim);
    VisItSetGetVariable(SimGetVariable, (void*)sim);
    VisItSetGetDomainNesting(SimGetDomainNesting, (void*)sim);
#ifdef PARALLEL
     VisItSetGetDomainList(SimGetDomainList, (void*)sim);
#endif
     
     while(!sim->done)
     {
        simulate_one_timestep(sim);
        VisItTimeStepChanged();

        // add a pseudocolor plot
        VisItAddPlot("Pseudocolor", "mandelbrot");
        VisItDrawPlots();
        VisItUpdatePlots();

        // save file
        char filename[100];
        sprintf(filename, "mandelbrot_amr_%06d.png", sim->saveCounter);
        if(VisItSaveWindow(filename, 1000, 1000, VISIT_IMAGEFORMAT_PNG) == VISIT_OKAY)
        {
            printf("Rendered to %s\n", filename);
        }
        else
        {
            printf("Failed to save %s\n", filename);
        }

        if( (sim->saveCounter % 10 == 0))
        {
            // export data to every 10th cycle
            visit_handle vars = VISIT_INVALID_HANDLE;
            VisIt_NameList_alloc(&vars);
            VisIt_NameList_addName(vars, "default");
            sprintf(filename, "mandelbrot_amr_export_%06d", sim->saveCounter);
            VisItExportDatabase(filename, "VTK_1.0", vars);
            printf("Exported to %s\n", filename);
            VisIt_NameList_free(vars);
        }
        sim->saveCounter++;
        VisItDeleteActivePlots();
     }
}

/******************************************************************************
 *
 * Purpose: This is the main function for the program.
 *
 * Programmer: Brad Whitlock
 * Date:       Thu Mar 19 11:54:04 PDT 2009
 *
 * Input Arguments:
 *   argc : The number of command line arguments.
 *   argv : The command line arguments.
 *
 * Modifications:
 *  Cyrus Harrison,Fri Feb  2 14:06:33 PST 2018
 *  Modify to change into batch example
 * 
 *****************************************************************************/

int main(int argc, char **argv)
{
    char *env = NULL;
    simulation_data sim;
    simulation_data_ctor(&sim);

#ifdef PARALLEL
    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &sim.par_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &sim.par_size);
#else
    sim.par_rank = 0;
    sim.par_size = 1;
#endif

    /* Initialize environment variables. */
    SimulationArguments(argc, argv);
    VisItSetupEnvironment();

    for(int i = 1; i < argc; ++i)
    {
        if(strcmp(argv[i], "-echo") == 0)
            sim.echo = true;
    }

#ifdef PARALLEL
    /* Install callback functions for global communication. */
    VisItSetBroadcastIntFunction2(visit_broadcast_int_callback, (void*)&sim);
    VisItSetBroadcastStringFunction2(visit_broadcast_string_callback, (void*)&sim);

    /* Tell libsim whether the simulation is parallel. */
    VisItSetParallel(sim.par_size > 1);
    VisItSetParallelRank(sim.par_rank);
#endif
    /* Only read the environment on rank 0. This could happen before MPI_Init if
     * we are using an MPI that does not like to let us spawn processes but we
     * would not know our processor rank.
     */
    if(sim.par_rank == 0)
        env = VisItGetEnvironment();

    /* Pass the environment to all other processors collectively. */
    VisItSetupEnvironment2(env);
    if(env != NULL)
        free(env);

    /* Init libsim. */
    if(sim.par_rank == 0)
    {
        VisItInitializeRuntime();
    }

    /* Read input problem setup, geometry, data. */
    read_input_deck();

    /* Call the main loop. */
    mainloop_batch(&sim);

    simulation_data_dtor(&sim);
#ifdef PARALLEL
    MPI_Finalize();
#endif

    return 0;
}

/* DATA ACCESS FUNCTIONS */

/******************************************************************************
 *
 * Purpose: This callback function returns simulation metadata.
 *
 * Programmer: Brad Whitlock
 * Date:       Thu Mar 19 11:54:04 PDT 2009
 *
 * Modifications:
 *
 *****************************************************************************/

visit_handle
SimGetMetaData(void *cbdata)
{
    visit_handle md = VISIT_INVALID_HANDLE;
    simulation_data *sim = (simulation_data *)cbdata;

    /* Create metadata. */
    if(VisIt_SimulationMetaData_alloc(&md) == VISIT_OKAY)
    {
        visit_handle mmd = VISIT_INVALID_HANDLE;
        visit_handle vmd = VISIT_INVALID_HANDLE;

        /* Set the simulation state. */
        VisIt_SimulationMetaData_setMode(md, (sim->runMode == SIM_STOPPED) ?
            VISIT_SIMMODE_STOPPED : VISIT_SIMMODE_RUNNING);
        VisIt_SimulationMetaData_setCycleTime(md, sim->cycle, sim->time);

        /* Fill in the AMR metadata. */
        if(VisIt_MeshMetaData_alloc(&mmd) == VISIT_OKAY)
        {
            /* Set the mesh's properties.*/
            VisIt_MeshMetaData_setName(mmd, "AMR_mesh");
            VisIt_MeshMetaData_setMeshType(mmd, VISIT_MESHTYPE_AMR);
            VisIt_MeshMetaData_setTopologicalDimension(mmd, 2);
            VisIt_MeshMetaData_setSpatialDimension(mmd, 2);

            int ndoms = patch_num_patches(&sim->patch);
            VisIt_MeshMetaData_setNumDomains(mmd, ndoms);
            VisIt_MeshMetaData_setDomainTitle(mmd, "patches");
            VisIt_MeshMetaData_setDomainPieceName(mmd, "patch");

            int nlevels = patch_num_levels(&sim->patch);
            VisIt_MeshMetaData_setNumGroups(mmd, nlevels);
            VisIt_MeshMetaData_setGroupTitle(mmd, "levels");
            VisIt_MeshMetaData_setGroupPieceName(mmd, "level");
            patch_t **patches = patch_flat_array(&sim->patch);
            int *pcount = (int *)malloc(nlevels * sizeof(int));
            memset(pcount, 0, nlevels * sizeof(int));
            for(int i = 0; i < ndoms; ++i)
            {
                char tmpName[100];
                sprintf(tmpName, "level%d,patch%04d", patches[i]->level, pcount[patches[i]->level]++);
                VisIt_MeshMetaData_addDomainName(mmd, tmpName);
                VisIt_MeshMetaData_addGroupId(mmd, patches[i]->level);
            }
            FREE(pcount);
            FREE(patches);

            VisIt_MeshMetaData_setXUnits(mmd, "cm");
            VisIt_MeshMetaData_setYUnits(mmd, "cm");
            VisIt_MeshMetaData_setXLabel(mmd, "Real");
            VisIt_MeshMetaData_setYLabel(mmd, "Imaginary");

            VisIt_SimulationMetaData_addMesh(md, mmd);
        }

        /* Add a variable. */
        if(VisIt_VariableMetaData_alloc(&vmd) == VISIT_OKAY)
        {
            VisIt_VariableMetaData_setName(vmd, "mandelbrot");
            VisIt_VariableMetaData_setMeshName(vmd, "AMR_mesh");
            VisIt_VariableMetaData_setType(vmd, VISIT_VARTYPE_SCALAR);
            VisIt_VariableMetaData_setCentering(vmd, VISIT_VARCENTERING_ZONE);

            VisIt_SimulationMetaData_addVariable(md, vmd);
        }

        /* Add some custom commands. */
        for(size_t i = 0; i < sizeof(cmd_names)/sizeof(const char *); ++i)
        {
            visit_handle cmd = VISIT_INVALID_HANDLE;
            if(VisIt_CommandMetaData_alloc(&cmd) == VISIT_OKAY)
            {
                VisIt_CommandMetaData_setName(cmd, cmd_names[i]);
                VisIt_SimulationMetaData_addGenericCommand(md, cmd);
            }
        }
    }

    return md;
}

/******************************************************************************
 *
 * Purpose: This callback function returns meshes.
 *
 * Programmer: Brad Whitlock
 * Date:       Thu Mar 19 11:54:04 PDT 2009
 *
 * Modifications:
 *
 *****************************************************************************/

visit_handle
SimGetMesh(int domain, const char *name, void *cbdata)
{
    visit_handle h = VISIT_INVALID_HANDLE;

    /* Get the patch with the appropriate domain id. */
    simulation_data *sim = (simulation_data *)cbdata;
    patch_t *patch = patch_get_patch(&sim->patch, domain);

    if(strcmp(name, "AMR_mesh") == 0 && patch != NULL)
    {
        int i;

        if(VisIt_RectilinearMesh_alloc(&h) != VISIT_ERROR)
        {
            /* Initialize X coords. */
            float *coordX = (float *)malloc(sizeof(float) * (patch->nx+1));
            float width  = sim->patch.window[1] - sim->patch.window[0];
            float x0 = (patch->window[0] - sim->patch.window[0]) / width;
            float x1 = (patch->window[1] - sim->patch.window[0]) / width;
            for(i = 0; i < (patch->nx+1); ++i)
            {
                float t = float(i) / float(patch->nx);
                coordX[i] = (1.f-t)*x0 + t*x1;
            }
            /* Initialize Y coords. */
            float *coordY = (float *)malloc(sizeof(float) * (patch->ny+1));
            float height = sim->patch.window[3] - sim->patch.window[2];
            float y0 = (patch->window[2] - sim->patch.window[2]) / height;
            float y1 = (patch->window[3] - sim->patch.window[2]) / height;
            for(i = 0; i < (patch->ny+1); ++i)
            {
                float t = float(i) / float(patch->ny);
                coordY[i] = (1.f-t)*y0 + t*y1;
            }

            /* Give the mesh some coordinates it can use. */
            visit_handle xc, yc;
            VisIt_VariableData_alloc(&xc);
            VisIt_VariableData_alloc(&yc);
            if(xc != VISIT_INVALID_HANDLE && yc != VISIT_INVALID_HANDLE)
            {
                VisIt_VariableData_setDataF(xc, VISIT_OWNER_VISIT, 1, patch->nx+1, coordX);
                VisIt_VariableData_setDataF(yc, VISIT_OWNER_VISIT, 1, patch->ny+1, coordY);
                VisIt_RectilinearMesh_setCoordsXY(h, xc, yc);
            }
            else
            {
                free(coordX);
                free(coordY);
            }
        }
    }

    return h;
}

/******************************************************************************
 *
 * Purpose: This callback function returns the domain nesting for the mesh.
 *
 * Programmer: Brad Whitlock
 * Date:       Thu Mar 19 14:31:51 PDT 2009
 *
 * Modifications:
 *   Brad Whitlock, Mon Aug  6 14:56:11 PDT 2012
 *   Added debug printing.
 *
 *****************************************************************************/

visit_handle
SimGetDomainNesting(const char *name, void *cbdata)
{
    simulation_data *sim = (simulation_data *)cbdata;
    visit_handle h = VISIT_INVALID_HANDLE;

    if(VisIt_DomainNesting_alloc(&h) != VISIT_ERROR)
    {
        int npatches = patch_num_patches(&sim->patch);
        int nlevels = patch_num_levels(&sim->patch);
        VisIt_DomainNesting_set_dimensions(h, npatches, nlevels, 2);
#ifdef DEBUG_PRINT
        std::ostringstream s;
        s << "Domain Nesting (cycle " << sim->cycle << "): npatches=" << npatches
          << ", nlevels=" << nlevels << ", ndims=2" << std::endl;
#endif

        int ratios[3];
        for(int level = 0; level < nlevels; ++level)
        {
            if(level == 0)
            {
                ratios[0] = ratios[1] = ratios[2] = 1;
            }
            else
            {
                ratios[0] = sim->refinement_ratio;
                ratios[1] = sim->refinement_ratio;
                ratios[2] = 1;
            }
            VisIt_DomainNesting_set_levelRefinement(h, level, ratios);
#ifdef DEBUG_PRINT
            s << "\tlevel " << level << " refinement " << ratios[0] << ", " << ratios[1] << ", " << ratios[2] << std::endl;
#endif
        }

        patch_t **patches = patch_flat_array(&sim->patch);
        int logicalExtents[6];
        for(int i = 0; i < npatches; ++i)
        {
            int *child_patches = NULL;
            int nchild_patches = patches[i]->nsubpatches;
            if(nchild_patches > 0)
            {
                child_patches = ALLOC(nchild_patches, int);
                for(int j = 0; j < nchild_patches; ++j)
                    child_patches[j] = patches[i]->subpatches[j].id;
            }
#ifdef DEBUG_PRINT
            s << "\tpatch " << i << " id=" << patches[i]->id << ", level=" << patches[i]->level << std::endl;
            s << "\t\txext=" << patches[i]->logical_extents[0] << ", " << patches[i]->logical_extents[1] << std::endl;
            s << "\t\tyext=" << patches[i]->logical_extents[2] << ", " << patches[i]->logical_extents[3] << std::endl;
            s << "\t\tnChildren=" << nchild_patches << ", children={";
            for(int j = 0; j < nchild_patches; ++j)
                s << child_patches[j] << ", ";
            s << "}" << std::endl;
#endif
            // logical extents are stored lowI,lowJ,lowL,hiI,hiJ,hiK
            logicalExtents[0] = patches[i]->logical_extents[0];
            logicalExtents[1] = patches[i]->logical_extents[2];
            logicalExtents[2] = 0;
            logicalExtents[3] = patches[i]->logical_extents[1];
            logicalExtents[4] = patches[i]->logical_extents[3];
            logicalExtents[5] = 0;
    
            VisIt_DomainNesting_set_nestingForPatch(h, 
                patches[i]->id, patches[i]->level, child_patches, nchild_patches,
                logicalExtents);
    
            FREE(child_patches);
        }
        FREE(patches);

#ifdef DEBUG_PRINT
        if(sim->par_rank == 0)
            std::cout << s.str();
#endif
    }

    return h;
}

/******************************************************************************
 *
 * Purpose: This callback function returns scalars.
 *
 * Programmer: Brad Whitlock
 * Date:       Thu Mar 19 14:31:51 PDT 2009
 *
 * Modifications:
 *
 *****************************************************************************/

visit_handle
SimGetVariable(int domain, const char *name, void *cbdata)
{
    visit_handle h = VISIT_INVALID_HANDLE;

    /* Get the patch with the appropriate domain id. */
    simulation_data *sim = (simulation_data *)cbdata;
    patch_t *patch = patch_get_patch(&sim->patch, domain);

    if(strcmp(name, "mandelbrot") == 0 && patch != NULL)
    { 
        VisIt_VariableData_alloc(&h);
        VisIt_VariableData_setDataC(h, VISIT_OWNER_SIM, 1,
            patch->nx * patch->ny, (char *)patch->data);
    }

    return h;
}

/******************************************************************************
 *
 * Purpose: This callback function returns a domain list.
 *
 * Programmer: Brad Whitlock
 * Date:       Mon Aug  6 14:55:17 PDT 2012
 *
 * Modifications:
 *
 *****************************************************************************/

visit_handle
SimGetDomainList(const char *name, void *cbdata)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    if(VisIt_DomainList_alloc(&h) != VISIT_ERROR)
    {
        visit_handle hdl = VISIT_INVALID_HANDLE;
        simulation_data *sim = (simulation_data *)cbdata;
        int npatches = patch_num_patches(&sim->patch);

        if(sim->npatch_list > 0)
        {
            VisIt_VariableData_alloc(&hdl);
            VisIt_VariableData_setDataI(hdl, VISIT_OWNER_SIM, 1, sim->npatch_list, sim->patch_list);
        }

        VisIt_DomainList_setDomains(h, npatches, hdl);
    }
    return h;
}
