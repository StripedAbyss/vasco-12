#pragma once
#ifndef __GENERALMESH_H__
#define __GENERALMESH_H__
#include <vector>
#include <fstream>
#include <string>
#include <Eigen/Dense>

using namespace std;
typedef vector<int> v_int;
const double PI = 3.14159265358979323;

class General_Mesh
{
public:
	General_Mesh();
	~General_Mesh();




	void genResultMesh(const char* filename);
	void genResultMesh_2(vector<int> num, vector<vector<double>> colors,const char* filename);
	void insert(vector<Eigen::Vector3d>& origin_vertexs, vector<v_int>& origin_faces,
		vector<Eigen::Vector3d> vertexs, vector<v_int> faces);
	void insert(vector<Eigen::Vector3d> vertexs, vector<v_int> faces);

	void genCube(double half_length);
	void genCylinder();

	void genHollowCylinder(float inside_radius, float outside_radius, float height);
	void genHollowCube(float width, float radius, float thick);



	vector<Eigen::Vector3d> meshTrans(Eigen::Vector3d dire, vector<Eigen::Vector3d> point);
	vector<Eigen::Vector3d> meshScale(Eigen::Vector3d scale, vector<Eigen::Vector3d> point);
	vector<Eigen::Vector3d> meshRotate(Eigen::Vector3d z_dire, vector<Eigen::Vector3d> point);
	vector<Eigen::Vector3d> meshRotate(Eigen::Vector3d x_dire, Eigen::Vector3d z_dire, vector<Eigen::Vector3d> point);



	vector<v_int> result_faces;
	vector<Eigen::Vector3d> result_vertex;
	vector<v_int> cylinder_faces;
	vector<Eigen::Vector3d> cylinder_vertex;
	vector<v_int> cube_faces;
	vector<Eigen::Vector3d> cube_vertex;
	vector<v_int> hollow_cylinder_faces;
	vector<Eigen::Vector3d> hollow_cylinder_vertex;
	vector<v_int> hollow_cube_faces;
	vector<Eigen::Vector3d> hollow_cube_vertex;
	double r;
	double g;
	double b;

};

#endif