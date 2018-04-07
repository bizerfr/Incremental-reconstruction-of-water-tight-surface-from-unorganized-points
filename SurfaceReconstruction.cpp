// SurfaceReconstruction.cpp : contain the main function
#include "stdafx.h"
#include<windows.h>
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkDataSetMapper.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkGlyph3D.h>
#include <vtkIdFilter.h>
#include <vtkLabeledDataMapper.h>
#include <vtkActor2D.h>
#include <vtkSTLWriter.h>
#include <vtksys/SystemTools.hxx>
#include <vtkObjectFactory.h>
#include <vtkProperty.h>
#include <vtkWindowToImageFilter.h>
#include <vtkJPEGWriter.h>
#include <vtkPNGWriter.h>
#include <vtkBMPWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkRenderWindowInteractor.h>
#include "SurfaceDelaunay.h"

using namespace std;


vtkRenderWindow *renWin = vtkRenderWindow::New();
//xyz file path
string pointCloudFile_path("data/pear_r.xyz");

// Define interaction style
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:

	static KeyPressInteractorStyle* New();
	vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);

	virtual void OnKeyPress() 
	{
		// Get the keypress
		vtkRenderWindowInteractor *rwi = this->Interactor;
		std::string key = rwi->GetKeySym();

		// Output the key that was pressed
		std::cout << "Pressed " << key << std::endl;


		// Handle a "normal" key
		if(key == "s")
		{
			std::cout << "screenshot" << std::endl;

			auto screenshot=[]()
			{
				// Screenshot
				vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = 
					vtkSmartPointer<vtkWindowToImageFilter>::New();
				windowToImageFilter->SetInput(renWin);
				//set the resolution of the output image (3 times the current resolution of vtk render window)
				windowToImageFilter->SetMagnification(5); 
				//also record the alpha (transparency) channel
				windowToImageFilter->SetInputBufferTypeToRGBA(); 
				// read from the back buffer
				windowToImageFilter->ReadFrontBufferOff(); 
				windowToImageFilter->Update();

				vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();

				int slash_pos = pointCloudFile_path.find_last_of('/');
				int dot_pos = pointCloudFile_path.find_last_of('.');
				string outputName = pointCloudFile_path.substr(slash_pos + 1, dot_pos - slash_pos -1 );
				string output_path=pointCloudFile_path.substr(0,4);
				output_path=output_path + "/" + outputName + ".tiff";

				writer->SetFileName(output_path.c_str());

				writer->SetInputConnection(windowToImageFilter->GetOutputPort());
				writer->Write();


			};

			screenshot();			

		}

		// Forward events
		vtkInteractorStyleTrackballCamera::OnKeyPress();
	}

};

vtkStandardNewMacro(KeyPressInteractorStyle);

int _tmain(int argc, _TCHAR* argv[])
{
	// true: insert isolated vertices; false: ignore isolated vertices;default:true;
	bool insert_iso = true;

	//-----can be commented out in DEBUG model-------
	if ( argc > 3 )
	{
		cout << "please input: SurfaceReconstruction xyz_file_path" 
			"T/F for inserting isolated vertices (default value is true)" << endl;
		exit(1);
	}
	else if (argc == 2)
	{
		pointCloudFile_path = argv[1];
	}
	else if (argc == 3)
	{
		pointCloudFile_path = argv[1];
		if (strcmp(argv[2],"T") == 0)
		{
			cout << "Insert isolated vertices" << endl;
			insert_iso = true;
		}
		else if (strcmp(argv[2], "F") ==0 )
		{
			cout << "Ignore isolated vertices" << endl;
			insert_iso = false;
		}
		else
		{
			cout << "please input T,F or defalut  for the 3rd argument" << endl;
			exit(1);
		}
	}
	//-----can be commented out in DEBUG model-------

	

    srand((unsigned)time(NULL));
	//system("mode con:cols=100 lines=10000");


	//store the raw points in the memory "pointsRaw"
	DataArray<double> pointsRaw(3,0);
	
	ifstream infile(pointCloudFile_path.c_str());//open xyz file
	if (!infile)
	{
		cout << "Unable to input nodes!";
		exit(1);
	}
	string line;

	int ii = 0;
	while (getline(infile, line))
	{
		auto first_string_space = line.find_first_of(' ');
		auto second_string_space = line.find_first_of(' ', first_string_space + 1 );
		double xx, yy, zz;
		if ((first_string_space == string::npos) | (second_string_space == string::npos))
		{
			cout << "incorrect .xyz file" << endl;
		}
		else
		{
			auto sieze_yy = second_string_space - first_string_space - 1;
			string tmp = line.substr(0, first_string_space );
			xx = atof(tmp.c_str());
			tmp = line.substr(first_string_space + 1, sieze_yy);
			yy = atof(tmp.c_str());
			tmp = line.substr(second_string_space + 1);
			zz = atof(tmp.c_str());
			double ptr[3];
			ptr[0] = xx;
			ptr[1] = yy;
			ptr[2] = zz;
			
			pointsRaw.insert_next_tuple(ptr);
			ii++;
		}
	}
	infile.close();

	//---------------bounding box----------------//
	double b0[3] = {-3000,-3000,-3000};
	double b1[3] = {3000,3000,11000};

	double **bounds;
	bounds=new double*[2];
	bounds[0] = b0;
	bounds[1] = b1;
	


	time_t t_start, t_end1,t_end2;
	t_start =  clock();

	//initialization for the constructor
	SurfaceDelaunay<double, double> Del(3,bounds,true);
	Del.insert_iso = insert_iso; // true: insert isolated vertices; false: ignore isolated vertices
	
	//insert the point by one one
	for(int i=0;i<=ii-1;i++)
	{
		const double* ptr = pointsRaw.get_tuple(i);

		//dynamic reconstruction
		Del.insert_point(ptr);
		//if (i==215)
		//{
		//	break;
		//}
	}
	//re-insert temporally abandoned points
	Del.handle_candidate_points();
	t_end1 = clock();

	cout<<endl<< "insert time:" <<double( t_end1 - t_start )/CLOCKS_PER_SEC<<" second"<< endl;



	DataArray<Idtype> surfaceData;
	DataArray<double> pointSet;
	surfaceData.init(3,0);
	pointSet.init(3,0);

	//output the vertices and surface facets for rendering
	Del.get_surface_data_structure(pointSet,surfaceData);


	// this part is used for rendering
	vtkSmartPointer<vtkPolyData> surface =
		vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	
	vtkSmartPointer<vtkCellArray> polys =
		vtkSmartPointer<vtkCellArray>::New();


	//get the points in the surface mesh
	for(int i=1;i<= pointSet.get_max_tuple_id();i++)
	{
		const double* ptrs= pointSet.get_tuple(i);
		points->InsertPoint(i,ptrs);    //no need to start from 0 in VTK
	}

	//get surface facets in the surface mesh
	for(int j=0;j<=surfaceData.get_max_tuple_id();j++)
	{
		vtkIdType aa[3] = { surfaceData.get_tuple_element(j,0), 
			 surfaceData.get_tuple_element(j,1),  
			 surfaceData.get_tuple_element(j,2) };
		const vtkIdType* fl = aa;
		polys->InsertNextCell(3, fl);
	}


	surface->SetPoints(points);
	points->Delete();
	surface->SetPolys(polys);
	polys->Delete();


	//write stl
	vtkSmartPointer<vtkSTLWriter> stlWriter =
		vtkSmartPointer<vtkSTLWriter>::New();
	stlWriter->SetInputData(surface);
	int slash_pos = pointCloudFile_path.find_last_of('/');
	int dot_pos = pointCloudFile_path.find_last_of('.');
	string outputStlName = pointCloudFile_path.substr(slash_pos + 1, dot_pos - slash_pos - 1);
	string outputStl_path=pointCloudFile_path.substr(0,4);
	outputStl_path=outputStl_path + "/" + outputStlName+".stl";
	stlWriter->SetFileName(outputStl_path.c_str());
	stlWriter->Write();


	vtkSmartPointer<vtkPolyDataMapper> surfaceMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
    surfaceMapper->SetInputData(surface);
	vtkSmartPointer<vtkActor> surfaceActor =
		vtkSmartPointer<vtkActor>::New();
	surfaceActor->GetProperty()->SetColor(1.0000, 0.8, 0);
    surfaceActor->SetMapper(surfaceMapper);
	vtkSmartPointer<vtkRenderer> ren =
		vtkSmartPointer<vtkRenderer>::New();
	renWin->AddRenderer(ren);
	renWin->SetSize(800,800);
	//interactor style
	vtkSmartPointer<vtkRenderWindowInteractor> iren =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);
	vtkSmartPointer<KeyPressInteractorStyle> keyPress =
		vtkSmartPointer<KeyPressInteractorStyle>::New();
	iren->SetInteractorStyle(keyPress);
	keyPress->SetCurrentRenderer(ren);
	ren->AddActor(surfaceActor);
	ren->SetBackground(1,1,1);
	renWin->Render();
	iren->Start();

	return EXIT_SUCCESS;
}

