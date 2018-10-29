// #include <openvdb/openvdb.h>
// #include <iostream>
// int main()
// {
//     // Initialize the OpenVDB library.  This must be called at least
//     // once per program and may safely be called multiple times.
//     openvdb::initialize();
//     // Create an empty floating-point grid with background value 0.
//     openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
//     std::cout << "Testing random access:" << std::endl;
//     // Get an accessor for coordinate-based access to voxels.
//     openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
//     // Define a coordinate with large signed indices.
//     openvdb::Coord xyz(1000, -200000000, 30000000);
//     // Set the voxel value at (1000, -200000000, 30000000) to 1.
//     accessor.setValue(xyz, 1.0);
//     // Verify that the voxel value at (1000, -200000000, 30000000) is 1.
//     std::cout << "Grid" << xyz << " = " << accessor.getValue(xyz) << std::endl;
//     // Reset the coordinates to those of a different voxel.
//     xyz.reset(1000, 200000000, -30000000);
//     // Verify that the voxel value at (1000, 200000000, -30000000) is
//     // the background value, 0.
//     std::cout << "Grid" << xyz << " = " << accessor.getValue(xyz) << std::endl;
//     // Set the voxel value at (1000, 200000000, -30000000) to 2.
//     accessor.setValue(xyz, 2.0);
//     // Set the voxels at the two extremes of the available coordinate space.
//     // For 32-bit signed coordinates these are (-2147483648, -2147483648, -2147483648)
//     // and (2147483647, 2147483647, 2147483647).
//     accessor.setValue(openvdb::Coord::min(), 3.0f);
//     accessor.setValue(openvdb::Coord::max(), 4.0f);
//     std::cout << "Testing sequential access:" << std::endl;
//     // Print all active ("on") voxels by means of an iterator.
//     for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter) {
//         std::cout << "Grid" << iter.getCoord() << " = " << *iter << std::endl;
//     }
// }

#include <iostream>
#include <fstream>
#include<sstream>
#include <vector>
#include <openvdb/openvdb.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointCount.h>
#include <string>

using namespace std;
int main()
{
    // Initialize grid types and point attributes types.
    openvdb::initialize();
    // Create a vector with four point positions.
    std::vector<openvdb::Vec3R> positions;
    int countT = 0;
    string filename = "/home/zhanghaitong/Desktop/vdbpoint.txt";
    

    std::ifstream inf(filename);
    if(!inf.is_open())
    {
        std::cout<< "error open file";
    }
    float number1,number2,number3;
    int cco = 0;
    while(!inf.eof())
    {
        
        //cout<<"in while"<<endl;
        inf >> number1 >> number2  >> number3 ;
      
        cout<<"number1: "<<number1<<cco<<endl;

        if(number1>96999.9)
        {
            cco+=1;
        }
        if((number1 < 99699.9) && (number2 < 99699.9) && (number3 < 99699.9) )
        {
           // cout<<
            //cout<<"in &&&&&&&&&&&&&&&&&"<<endl;
            positions.push_back(openvdb::Vec3R(number1, number2, number3));
        }
        else
        {

            countT++;
            cout << countT << endl;
            openvdb::points::PointAttributeVector<openvdb::Vec3R> positionsWrapper(positions);
            int pointsPerVoxel = 8;
            float voxelSize =
                openvdb::points::computeVoxelSize(positionsWrapper, pointsPerVoxel);
            // Print the voxel-size to cout
            std::cout << "VoxelSize=" << voxelSize << std::endl;
            // Create a transform using this voxel-size.
            openvdb::math::Transform::Ptr transform =
                openvdb::math::Transform::createLinearTransform(voxelSize);
            // Create a PointDataGrid containing these four points and using the
            // transform given. This function has two template parameters, (1) the codec
            // to use for storing the position, (2) the grid we want to create
            // (ie a PointDataGrid).
            // We use no compression here for the positions.
            openvdb::points::PointDataGrid::Ptr grid =
                openvdb::points::createPointDataGrid<openvdb::points::NullCodec,
                            openvdb::points::PointDataGrid>(positions, *transform);
            // Set the name of the grid
            grid->setName("Points");
            // Create a VDB file object for writing.
            string name = to_string(countT);

            string pointname = "/home/zhanghaitong/Desktop/ttt/points_data/mypoint";
            int k = 5-name.size();
            while(k--)
            {
                name = "0"+ name;
            }
            pointname += name;
            pointname += ".vdb";
            std::cout<<name<<"zzzz"<<std::endl;
            
            openvdb::io::File file(pointname);
            // Add the PointDataGrid pointer to a container.
            openvdb::GridPtrVec grids;
            grids.push_back(grid);
            // Write out the contents of the container.
            file.write(grids);
            file.close();
            // Create a new VDB file object for reading.
            openvdb::io::File newFile(pointname);
            // Open the file. This reads the file header, but not any grids.
            newFile.open();
            // Read the grid by name.
            openvdb::GridBase::Ptr baseGrid = newFile.readGrid("Points");
            newFile.close();
            // From the example above, "Points" is known to be a PointDataGrid,
            // so cast the generic grid pointer to a PointDataGrid pointer.
            grid = openvdb::gridPtrCast<openvdb::points::PointDataGrid>(baseGrid);
            openvdb::Index64 count = openvdb::points::pointCount(grid->tree());
            std::cout << "PointCount=" << count << std::endl;
            positions.clear();
        }
    }
}