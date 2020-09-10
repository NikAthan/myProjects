#include "PointClouds.h"
#include<iostream>
#include <windows.h>
#include <stdio.h> 
#include <string>
#include <limits>
#include <Eigen\Dense>
#include <opencv2/imgproc/imgproc.hpp>
#include <chrono>
#define DIMENSIONS 3
#undef max

int l0,l1,l2,l3;
int q0 = 0; int q1 = 1; int q0ii = 1;
int N;int Ni;int Nj; int N5;
bool er1 = true; bool er2 = false;bool er4 = false;
bool er5 = false; bool er6 = false; bool er7 = false;
bool er8 = false; bool er9 = false;
int D = 10;  //  Βήμα Δειγματοληψίας
 

#define SPLIT_INSTEAD_OF_INTERSECT 0

using namespace std;
using namespace vvr;
using namespace Eigen;



//! KDTree::

KDTree::KDTree(VecVector &pts)
	: pts(pts)
{
	const float t = vvr::getSeconds();
	m_root = new KDNode();
	m_depth = makeNode(m_root, pts, 0);
	const float KDTree_construction_time = vvr::getSeconds() - t;
}

KDTree::~KDTree()
{
	const float t = vvr::getSeconds();
	delete m_root;
	const float KDTree_destruction_time = vvr::getSeconds() - t;
}

int KDTree::makeNode(KDNode *node, VecVector &pts, const int level)
{
	//! Sort along the appropriate axis, find median point and split.
	const int axis = level % DIMENSIONS;
	std::sort(pts.begin(), pts.end(), VecComparator(axis));
	const int i_median = pts.size() / 2;

	//! Set node members
	node->level = level;
	node->axis = axis;
	node->split_point = pts[i_median];
	node->aabb.SetFrom(&pts[0], pts.size());

	//! Continue recursively or stop.
	if (pts.size() <= 1)
	{
		return level;
	}
	else
	{
		int level_left = 0;
		int level_right = 0;

		VecVector pts_left(pts.begin(), pts.begin() + i_median);
		VecVector pts_right(pts.begin() + i_median + 1, pts.end());

		if (!pts_left.empty())
		{
			node->child_left = new KDNode();
			level_left = makeNode(node->child_left, pts_left, level + 1);

		}
		if (!pts_right.empty())
		{
			node->child_right = new KDNode();
			level_right = makeNode(node->child_right, pts_right, level + 1);
		}

		int max_level = std::max(level_left, level_right);
		return max_level;
	}
}

void KDTree::getNodesOfLevel(KDNode *node, std::vector<KDNode*> &nodes, int level)
{
	if (!level)
	{
		nodes.push_back(node);
	}
	else
	{
		if (node->child_left) getNodesOfLevel(node->child_left, nodes, level - 1);
		if (node->child_right) getNodesOfLevel(node->child_right, nodes, level - 1);
	}
}

std::vector<KDNode*> KDTree::getNodesOfLevel(const int level)
{
	std::vector<KDNode*> nodes;
	if (!m_root) return nodes;
	getNodesOfLevel(m_root, nodes, level);
	return nodes;
}

// Struct for Rotation-Translation

struct RotTran {
	MatrixXd R;
	Vector3d t;
};



RotTran rotran(vvr::PointCloud pointc1, vvr::PointCloud pointc2) {
	
	RotTran rtrn;
	double w = 0.1;
	int n = pointc1.size();
	PointCloud XW;
	PointCloud Y;

	//// STEP 1
	
	Point3D num1(0, 0, 0);
	Point3D num2(0, 0, 0);
	for (int i = 0; i < n; i++) {
		num1.x = num1.x + w*pointc1[i].x;
		num1.y = num1.y + w*pointc1[i].y;
		num1.z = num1.z + w*pointc1[i].z;
		num2.x = num2.x + w*pointc2[i].x;
		num2.y = num2.y + w*pointc2[i].y;
		num2.z = num2.z + w*pointc2[i].z;
	}
	
	//Ypologismos p mpara = P
	double den1 = n*w;
	Point3D P(num1.x / den1 , num1.y / den1 , num1.z / den1);

	//Ypologismos q mpara = Q
	double den2 = n*w;
	Point3D Q(num2.x / den2, num2.y / den2, num2.z / den2);
	
	
	//// STEPS 2-3
	
	for (int i = 0; i < n; i++) {

		// kataskeui tou pinaka X*W
		Point3D vx;
		vx.x = (pointc1[i].x - P.x)*w;
		vx.y = (pointc1[i].y - P.y)*w;
		vx.z = (pointc1[i].z - P.z)*w;
		XW.push_back(vx);

		// kataskeui tou pinaka Y
		Point3D vy;
		vy.x = pointc2[i].x - Q.x;
		vy.y = pointc2[i].y - Q.y;
		vy.z = pointc2[i].z - Q.z;
		Y.push_back(vy);
	}
	
	
	// Ypologismos pinaka S=X*W*Y.transpose() se Eigen

	Eigen::MatrixXd S(3, 3);
	double sum11=0; double sum21=0; double sum31=0;
	double sum12=0; double sum22=0; double sum32=0;
	double sum13=0; double sum23=0; double sum33=0;
	for (int i = 0; i < n; i++) {
		sum11 = sum11 + XW[i].x*Y[i].x;
		sum21 = sum21 + XW[i].y*Y[i].x;
		sum31 = sum31 + XW[i].z*Y[i].x;
		sum12 = sum12 + XW[i].x*Y[i].y;
		sum22 = sum22 + XW[i].y*Y[i].y;
		sum32 = sum32 + XW[i].z*Y[i].y;
		sum13 = sum13 + XW[i].x*Y[i].z;
		sum23 = sum23 + XW[i].y*Y[i].z;
		sum33 = sum33 + XW[i].z*Y[i].z;
	}


	S.row(0) << sum11, sum12, sum13;
	S.row(1) << sum21, sum22, sum23;
	S.row(2) << sum31, sum32, sum33;
	

	//// STEP 4

	// Eksagwgi pinakwn U kai V

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
	
	// Ypologismos pinaka Σ

	MatrixXd SS;
	SS.setZero(3,3);
	SS(0, 0) = 1;
	SS(1, 1) = 1;
	SS(2, 2) = (svd.matrixV()*(svd.matrixU().transpose())).determinant();

	// Ypologismos Rotate Re se Eigen::MatrixXd

	Eigen::MatrixXd Re = svd.matrixV() * SS * (svd.matrixU()).transpose();    
	//cout << "Rotation: " << Re << endl;

	// Metatropi P,Q apo Point3D se Eigen::Vector3d -> PE,QE

	Eigen::Vector3d PE(P.x, P.y, P.z);
	Eigen::Vector3d QE(Q.x, Q.y, Q.z);

	// Ypologismos translation te se Eigen::Vector3f

	Eigen::VectorXd te = QE- Re*PE;   
	//cout << "Translation: " << te << endl;

	// Return

	rtrn.R = Re;
	rtrn.t = te;
	
	return rtrn;;
}

PointCloudScene::PointCloudScene()
{
	//! Load settings.
	vvr::Shape::DEF_LINE_WIDTH = 4;
	vvr::Shape::DEF_POINT_SIZE = 4;
	m_perspective_proj = true;
	m_bg_col = Colour("768E77");
	m_point_color = Colour("454545");
	m_hide_sliders = false;
	m_hide_log = false;   //true

	// Load Object
	const string objDir = getBasePath() + "resources/obj/";
	const string objFile = objDir + "suzanne.obj";
	m_model_original = vvr::Mesh(objFile);

	// Load Images
	//const string imgDir = getBasePath() + "resources/images/";
	//const string colorImgFile = imgDir + "sample.png";
	//const string depthImgFile = imgDir + "sample_depth.png";
	//m_color_image = cv::imread(colorImgFile, CV_LOAD_IMAGE_ANYCOLOR);
	//m_depth_image = cv::imread(depthImgFile, CV_LOAD_IMAGE_ANYDEPTH);

	// Load RGB-D Images
	const string imgDirection = getBasePath() + "resources/rgbd-scenes/meeting_small/meeting_small_1/";

	// Load depth images
	float maxx = 0;
	for (int i = 0; i < 179; i++) {
		double maxDepth = 0;
		string name = "meeting_small_1_";
		string name2 = to_string(i + 1);
		string name3 = "_depth.png";
		string Name = name + name2 + name3;
		string depthim = imgDirection + Name;
		m_depth_image_new = cv::imread(depthim, CV_LOAD_IMAGE_ANYDEPTH);
		cv::minMaxLoc(m_depth_image_new, (double*)0, &maxDepth);   // maxDepth tis mias eikonas          
		if (maxx < maxDepth) {
			maxx = maxDepth;   // maxDepth olwn twn eikonwn
		}
		FI.push_back(m_depth_image_new);  // depth images
	}

	// Load rgb images
	for (int i = 0; i < 179; i++) {
		cv::Mat m_depth_image_norm;    //arxikopoiisi-ananewsi
		string name = "meeting_small_1_";
		string name2 = to_string(i + 1);
		string name3 = ".png";
		string Name = name + name2 + name3;
		string colorImgFile = imgDirection + Name;
		// normalization
		FI[i].convertTo(m_depth_image_norm, CV_8U, 255.0f / maxx);
		F.push_back(m_depth_image_norm);  // normalized depth images
		I.push_back(cv::imread(colorImgFile, CV_LOAD_IMAGE_ANYCOLOR));  // colour images
	}


	for (int i = 0; i < 50; i++) {          // GIA OLES TIS EIKONES <179 
		vvr::PointCloud m_model_points_new;
		m_depth_image_norm = F[i];
		std::vector<vec> v;
		PC(m_depth_image_norm, v);
		vall.push_back(v);                   // To vall exei ola ta vertices olwn twn eikonwn
		for (auto& pt : v)
		{
			m_model_points_new.push_back(Point3D(pt.x, pt.y, pt.z, vvr::Colour::black));
		}
		pointclouds.push_back(m_model_points_new);
	}


    reset();
}

void PointCloudScene::reset()
{
    Scene::reset();

    //! Define plane
    m_plane_d = 0;
    m_plane = Plane(vec(0, 1, 1).Normalized(), m_plane_d);

    //! Define what will be vissible by default
    m_style_flag = 0;
    //m_style_flag |= FLAG_SHOW_SOLID;
    m_style_flag |= FLAG_SHOW_WIRE;
    m_style_flag |= FLAG_SHOW_AXES;
    m_style_flag |= FLAG_SHOW_POINTS;
    //m_style_flag |= FLAG_SHOW_PLANE
}

void PointCloudScene::resize()
{
    //! By Making `first_pass` static and initializing it to true,
    //! we make sure that the if block will be executed only once.
    static bool first_pass = true;

    if (first_pass)
    {
        m_model_original.setBigSize(getSceneWidth() / 2);
        m_model_original.update();
        m_model = m_model_original;

		// Also get Point cloud from re-scaled model
		getPointCloud(m_model, m_model_points, vvr::Colour::green);

		// Initialize OpenCV Windows for the images
		cv::namedWindow("Color Image");
		cv::namedWindow("Depth Image");
		cv::namedWindow("Depthsssss Image");
		cv::namedWindow("Sobel");
		cv::namedWindow("Threshold");
        first_pass = false;
    }
}


void PointCloudScene::arrowEvent(ArrowDir dir, int modif)
{
    math::vec n = m_plane.normal;
    if (dir == UP) m_plane_d += 1;
    if (dir == DOWN) m_plane_d -= 1;
    else if (dir == LEFT) n = math::float3x3::RotateY(DegToRad(1)).Transform(n);
    else if (dir == RIGHT) n = math::float3x3::RotateY(DegToRad(-1)).Transform(n);
    m_plane = Plane(n.Normalized(), m_plane_d);

}

void PointCloudScene::keyEvent(unsigned char key, bool up, int modif)
{
    Scene::keyEvent(key, up, modif);
    key = tolower(key);

    switch (key)
    {
    case 's': m_style_flag ^= FLAG_SHOW_SOLID; break;
    case 'w': m_style_flag ^= FLAG_SHOW_WIRE; break;
    case 'n': m_style_flag ^= FLAG_SHOW_NORMALS; break;
    case 'a': m_style_flag ^= FLAG_SHOW_AXES; break;
    case 'p': m_style_flag ^= FLAG_SHOW_PLANE; break;
    case 'b': m_style_flag ^= FLAG_SHOW_POINTS; break;
	case 'd': (er2 = true); break;            // ERWTIMA DIO 2
	case 't': (er4 = true); break;            // ERWTIMA TESSERA 4
	case 'f': (er5 = true); break;            // ERWTIMA PENTE 5
	case 'e': (er6 = true); break;            // ERWTIMA EKSI 6
	case 'x': (er7 = true); break;            // ERWTIMA EFTA 7
	case 'l': (er8 = true); break;            // ERWTIMA PENTE 5 A
	case 'q': (er9 = true); break;            // ERWTIMA PENTE 5 B
    }
}

void PointCloudScene::draw()
{

	//! Draw plane
	if (m_style_flag & FLAG_SHOW_PLANE) {
		vvr::Colour colPlane(0x41, 0x14, 0xB3);
		float u = 20, v = 20;
		math::vec p0(m_plane.Point(-u, -v, math::vec(0, 0, 0)));
		math::vec p1(m_plane.Point(-u, v, math::vec(0, 0, 0)));
		math::vec p2(m_plane.Point(u, -v, math::vec(0, 0, 0)));
		math::vec p3(m_plane.Point(u, v, math::vec(0, 0, 0)));
		math2vvr(math::Triangle(p0, p1, p2), colPlane).draw();
		math2vvr(math::Triangle(p2, p1, p3), colPlane).draw();
	}
	/*
	if (m_style_flag & FLAG_SHOW_SOLID) m_model.draw(m_point_color, SOLID);
	if (m_style_flag & FLAG_SHOW_WIRE) m_model.draw(Colour::black, WIRE);
	if (m_style_flag & FLAG_SHOW_NORMALS) m_model.draw(Colour::black, NORMALS);
	if (m_style_flag & FLAG_SHOW_AXES) m_model.draw(Colour::black, AXES);
	if (m_style_flag & FLAG_SHOW_POINTS) drawPointCloud(m_model_points);
	*/


	if ((er2 == true) || (er4 == true) || (er5 == true) || (er6 == true) || (er7 == true)) {   // Na min ekteleitai to erwtima 1 otan ektelountai ta alla
		er1 = false;
	}
	if ((er1 == true) || (er4 == true) || (er5 == true) || (er6 == true) || (er7 == true)) {   // Na min ekteleitai to erwtima 2 otan ektelountai ta alla
		er2 = false;
		int f2 = 0;
	}
	if ((er1 == true) || (er2 == true) || (er5 == true) || (er6 == true) || (er7 == true)) {   // Na min ekteleitai to erwtima 4 otan ektelountai ta alla
		er4 = false;
		int f4 = 0;
	}
	if ((er1 == true) || (er2 == true) || (er4 == true) || (er6 == true) || (er7 == true)) {   // Na min ekteleitai to erwtima 5 otan ektelountai ta alla
		er5 = false;
		er8 = false;
		er9 = false;
		int f5 = 0;
	}
	if ((er1 == true) || (er2 == true) || (er4 == true) || (er5 == true) || (er7 == true)) {   // Na min ekteleitai to erwtima 6 otan ektelountai ta alla
		er6 = false;
		int f6 = 0;
	}
	if ((er1 == true) || (er2 == true) || (er4 == true) || (er5 == true) || (er6 == true)) {   // Na min ekteleitai to erwtima 7 otan ektelountai ta alla
		er7 = false;
		int f7 = 0;
	}
	if ((er8 == true) && (q0ii == 1)) {
		er9 = false;
	}
	if ((er9 == true) && (q1 == 1)) {
		er8 = false;
	}

	m_depth_image_norm = F[N];
	if ((er2 == true) || (er5 == true)) {
		if (N != 0) {
			m_depth_image_norm_previous = F[N - 1];
		}
	}
	m_color_image = I[N];
	m_depth = FI[N];

	///////////       ERWTIMA 2         

	if ((er2 == true) && (q0 == 1) && (N > 0)) {     // H sinthiki q0==1 gia na ekteleitai tin prwti fora kai mono otan allazw to N
		ls1.clear();
		// Plisiesteros geitonas gia kathe simeio
		int N1 = vall[N].size();
		int N2 = vall[N-1].size();
		float SinApostasi = 0;
		Point3D p;
		auto t0 = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < N1; i++) {      
			float minDistance = std::numeric_limits<float>::max();  // Arxikopoiisi tis apostasis ws apeiris
			for (int j = 0; j < N2; j++) {              
				float pointDistance = vall[N][i].DistanceSq(vall[N-1][j]);
				if (pointDistance < minDistance) {
					minDistance = pointDistance;
					p.x = vall[N - 1][j].x;
					p.y = vall[N - 1][j].y;
					p.z = vall[N - 1][j].z;
				}
			}
			SinApostasi = SinApostasi + minDistance;
			LineSeg3D lt(p.x, p.y, p.z, vall[N][i].x, vall[N][i].y, vall[N][i].z, vvr::Colour::blue);
			ls1.push_back(lt);
		}
		auto t1 = std::chrono::high_resolution_clock::now();
		auto dt = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
		cout << "XRONOS IPOLOGISMOY ANTISTOIXISIS SIMEIWN: " << dt << " seconds" << endl;
		float MesiApostasi = SinApostasi / N1;
		cout << "SINOLIKI APOSTASI: " << SinApostasi << endl;
		cout << "MESI APOSTASI: " << MesiApostasi << endl;
		q0 = 0;
	}

	
	///////////     ERWTIMA 4   

	if ((er4 == true) && (q1 == 1)) {
		int f = 0;
		float e = 0;
		vrn.clear();
		vtr.clear();
		auto t2 = std::chrono::high_resolution_clock::now();
		for (int i = Ni; i < Nj; i++) {
			float edistance = 0;
			if (f == 0) {     //   Na ekteleitai mono tin prwti fora - First Pass
				for (auto& v : vall[Ni])
				{
					vtr.push_back(v);
				}
				/// Euresi plisiesterwn vertices tou PC Ni ws pros to PC Ni+1
				std::vector<vec> P;
				float sinap = 0;
				float mesi = 0;
				Plisiesteros(vall[i + 1], vall[i], P, sinap, mesi);// P:plisiestera vertices tou PC i ws pros to PC i+1
				/// Eksagwgi PC plisiesterwn vertices
				vvr::PointCloud PNearest;
				getPointCloud2(P, PNearest);
				/// Rotation - Translation - Μετατόπιση και Περιστροφή καθενός Pointcloud(i+1) προς τα πλησιέστερα που ανήκουν στο Pointcloud(i)
				RotTran RT;
				RT = rotran(pointclouds[i + 1], PNearest);   // Μετατόπιση και Περιστροφή καθενός Pointcloud ως προς
				int M = pointclouds[i + 1].size();                      //   με τα πλησιέστερα του προηγούμενου
				for (int j = 0; j < M; j++) {
					float vv1 = (RT.R(0, 0)*pointclouds[i + 1][j].x + RT.R(0, 1)*pointclouds[i + 1][j].y + RT.R(0, 2)*pointclouds[i + 1][j].z + RT.t(0));
					float vv2 = (RT.R(1, 0)*pointclouds[i + 1][j].x + RT.R(1, 1)*pointclouds[i + 1][j].y + RT.R(1, 2)*pointclouds[i + 1][j].z + RT.t(1));
					float vv3 = (RT.R(2, 0)*pointclouds[i + 1][j].x + RT.R(2, 1)*pointclouds[i + 1][j].y + RT.R(2, 2)*pointclouds[i + 1][j].z + RT.t(2));
					Point3D pe(vv1, vv2, vv3);
					/// Αποστάσεις των ευθυγραμμισμένων από τα προηγούμενα(i) vertices
					edistance = edistance + (vall[i + 1][j].DistanceSq(P[j]));
					vrn.push_back(vec(vv1, vv2, vv3));  // Metatopismena vertices
				}
				/// Ως ανοχή e θεωρείται η μέση απόσταση των ευθυγραμμισμένων από τα προηγούμενα(i) vertices
				e = edistance / float(M);
				cout << "Anoxi ε: " << e << endl;
				/// Euresi plisiesterwn
				/// Εισαγωγή στο vtr των ευθυγραμμισμένων σημείων του i+1 ως προς το i που δεν είναι διπλά , δηλαδή σημείων που d > ε
				std::vector<vec> PP;
				int M1 = vtr.size();
				int M2 = vrn.size();
				for (int i = 0; i < M2; i++) {
					float d = std::numeric_limits<float>::max();
					vec n;
					for (int j = 0; j < M1; j++) {
						float Distancee = vtr[j].DistanceSq(vrn[i]);
						if (Distancee < d) {
							d = Distancee;
							n = vtr[j];
						}
					}
					float Distance = n.DistanceSq(vrn[i]);
					if (Distance > e) {
						vtr.push_back(vrn[i]);
					}
				}
			}
			else {
				/// Euresi plisiesterwn vertices tou PC i ws pros to PC i+1
				std::vector<vec> P;
				float sinap = 0;
				float mesi = 0;
				Plisiesteros(vall[i + 1], vrn, P, sinap, mesi);// P:plisiestera vertices tou PC i ws pros to PC i+1
				/// Eksagwgi PC plisiesterwn vertices
				vvr::PointCloud PNearest;
				getPointCloud2(P, PNearest);
				/// Rotation - Translation - Μετατόπιση και Περιστροφή καθενός Pointcloud(i+1) προς τα πλησιέστερα που ανήκουν στο Pointcloud(i)
				RotTran RT;
				RT = rotran(pointclouds[i + 1], PNearest);
				int M = pointclouds[i + 1].size();
				vrn.clear();
				for (int j = 0; j < M; j++) {
					float vv1 = (RT.R(0, 0)*pointclouds[i + 1][j].x + RT.R(0, 1)*pointclouds[i + 1][j].y + RT.R(0, 2)*pointclouds[i + 1][j].z + RT.t(0));
					float vv2 = (RT.R(1, 0)*pointclouds[i + 1][j].x + RT.R(1, 1)*pointclouds[i + 1][j].y + RT.R(1, 2)*pointclouds[i + 1][j].z + RT.t(1));
					float vv3 = (RT.R(2, 0)*pointclouds[i + 1][j].x + RT.R(2, 1)*pointclouds[i + 1][j].y + RT.R(2, 2)*pointclouds[i + 1][j].z + RT.t(2));
					Point3D pe(vv1, vv2, vv3);
					/// Εισαγωγή στο vrn των μετατοπισμένων σημείων
					vrn.push_back(vec(vv1, vv2, vv3));
				}
				/// Euresi plisiesterwn
				/// Εισαγωγή στο vtr των ευθυγραμμισμένων σημείων του i+1 ως προς το i που δεν είναι διπλά , δηλαδή σημείων που d > ε
				std::vector<vec> PP;
				int M1 = vtr.size();
				int M2 = vrn.size();
				for (int i = 0; i < M2; i++) {
					float d = std::numeric_limits<float>::max();
					vec n;
					for (int j = 0; j < M1; j++) {
						float Distancee = vtr[j].DistanceSq(vrn[i]);
						if (Distancee < d) {
							d = Distancee;
							n = vtr[j];
						}
					}
					float Distance = n.DistanceSq(vrn[i]);
					if (Distance > e) {
						vtr.push_back(vrn[i]);
					}
				}
				/// Αφαίρεση διπλών σημείων , δηλαδή σημείων που d < ε
				for (int j = 0; j < M2; j++) {
					
				}
			}
			f = 1;
		}
		cout << "Plithos vertices: " << vtr.size() << endl;
		getPointCloud2(vtr, pcM);
		auto t3 = std::chrono::high_resolution_clock::now();
		auto dt2 = std::chrono::duration_cast<std::chrono::seconds>(t3 - t2).count();
		cout << "XRONOS IPOLOGISMOY OLIKOY NEFOYS M: " << dt2 << " seconds" << endl;
		q1 = 0;
	}




	///////////    ERWTIMA 5


	///   Erwtima Aii me KDTREES

	if ((er5 == true) && (er8 == true) && (q0ii == 1) && (N5 > 0)) {
		ls.clear();
		auto t8 = std::chrono::high_resolution_clock::now();
		m_KDTree = new KDTree(vall[N5 - 1]);
		float best_dist;
		int S = vall[N5].size();
		for (int i = 0; i < S; i++) {
			best_dist = std::numeric_limits<float>::max();
			vec nn;
			const KDNode *nearest = NULL;
			NearestNeighboorKDtree(vall[N5][i], m_KDTree->root(), nearest, best_dist);
			if (nearest) nn = nearest->split_point; 
			LineSeg3D lt(nn.x, nn.y, nn.z, vall[N5][i].x, vall[N5][i].y, vall[N5][i].z, vvr::Colour::blue);
			ls.push_back(lt);
		}
		auto t9 = std::chrono::high_resolution_clock::now();
		auto dt5 = std::chrono::duration_cast<std::chrono::seconds>(t9 - t8).count();
		cout << "XRONOS Aii ME KDTREES: " << dt5 << " seconds" << endl;
		delete m_KDTree;
		q0ii = 0;
	}


	///   Erwtimata 3,4 me KDTREES

	if ((er5 == true) && (er9 == true) && (q1 == 1)){ 
		int f = 0;
		float e = 0;
		vrn.clear();
		vtr.clear();
		auto t4 = std::chrono::high_resolution_clock::now();
		for (int i = Ni; i < Nj; i++) {
			float edistance = 0;
			std::vector<vec> P;
			if (f == 0) {     //   Na ekteleitai mono tin prwti fora - First Pass
				for (auto& v : vall[Ni])
				{
					vtr.push_back(v);
				}
				/// Euresi plisiesterwn vertices tou PC Ni ws pros to PC Ni+1
				PlisiesterosKDTREE(vall[i + 1], vall[i], P);
				/// Eksagwgi PC plisiesterwn vertices
				vvr::PointCloud PNearest;
				getPointCloud2(P, PNearest);
				/// Rotation - Translation - Μετατόπιση και Περιστροφή καθενός Pointcloud(i+1) προς τα πλησιέστερα που ανήκουν στο Pointcloud(i)
				RotTran RT;
				RT = rotran(pointclouds[i + 1], PNearest);   // Μετατόπιση και Περιστροφή καθενός Pointcloud ως προς
				int M = pointclouds[i + 1].size();                      //   με τα πλησιέστερα του προηγούμενου
				for (int j = 0; j < M; j++) {
					float vv1 = (RT.R(0, 0)*pointclouds[i + 1][j].x + RT.R(0, 1)*pointclouds[i + 1][j].y + RT.R(0, 2)*pointclouds[i + 1][j].z + RT.t(0));
					float vv2 = (RT.R(1, 0)*pointclouds[i + 1][j].x + RT.R(1, 1)*pointclouds[i + 1][j].y + RT.R(1, 2)*pointclouds[i + 1][j].z + RT.t(1));
					float vv3 = (RT.R(2, 0)*pointclouds[i + 1][j].x + RT.R(2, 1)*pointclouds[i + 1][j].y + RT.R(2, 2)*pointclouds[i + 1][j].z + RT.t(2));
					Point3D pe(vv1, vv2, vv3);
					/// Αποστάσεις των πλησιέστερων από τα vertices του Ni+1 για υπολογισμό της ανοχής
					edistance = edistance + (vall[i+1][j].DistanceSq(P[j]));
					/// Εισαγωγή στο vrn των μετατοπισμένων σημείων
					vrn.push_back(vec(vv1, vv2, vv3));   //  Metatopismena vertices
				}
				/// Ως ανοχή e θεωρείται η μέση απόσταση των πλησιέστερων vertices tou PC Ni ws pros to PC Ni+1 από τα vertices του Ni+1
				/// Προσοχή!! Υπολογίζεται μόνο στο first Pass και στις επόμενες επαναλήψεις παραμένει το ίδιο 
				e = edistance / float(M);
				cout << "Anoxi ε: " << e << endl;
				/// Εισαγωγή στο vtr των ευθυγραμμισμένων σημείων του Νi+1 ως προς το Νi που δεν είναι διπλά , δηλαδή σημείων που d > ε
				KDTREEProsthesimidiplwn(vrn, vtr, e);
			}
			else {
				/// Euresi plisiesterwn vertices tou euthigrammismenou PC tis proigoumenis epanalipsis i ws pros to PC i+1
				PlisiesterosKDTREE(vall[i + 1], vrn, P);
				/// Eksagwgi PC plisiesterwn vertices
				vvr::PointCloud PNearest;
				getPointCloud2(P, PNearest);
				/// Rotation - Translation - Μετατόπιση και Περιστροφή καθενός Pointcloud(i+1) προς τα πλησιέστερα που ανήκουν στο Pointcloud(i)
				RotTran RT;
				RT = rotran(pointclouds[i + 1], PNearest);   // Μετατόπιση και Περιστροφή καθενός Pointcloud ως προς
				int M = pointclouds[i + 1].size();                      //   με τα πλησιέστερα του προηγούμενου ευθυγραμμισμένου
				vrn.clear();
				for (int j = 0; j < M; j++) {
					float vv1 = (RT.R(0, 0)*pointclouds[i + 1][j].x + RT.R(0, 1)*pointclouds[i + 1][j].y + RT.R(0, 2)*pointclouds[i + 1][j].z + RT.t(0));
					float vv2 = (RT.R(1, 0)*pointclouds[i + 1][j].x + RT.R(1, 1)*pointclouds[i + 1][j].y + RT.R(1, 2)*pointclouds[i + 1][j].z + RT.t(1));
					float vv3 = (RT.R(2, 0)*pointclouds[i + 1][j].x + RT.R(2, 1)*pointclouds[i + 1][j].y + RT.R(2, 2)*pointclouds[i + 1][j].z + RT.t(2));
					Point3D pe(vv1, vv2, vv3);
					/// Εισαγωγή στο vrn των μετατοπισμένων σημείων
					vrn.push_back(vec(vv1, vv2, vv3));
				}
				/// Εισαγωγή στο vtr των ευθυγραμμισμένων σημείων του i+1 ως προς το i που δεν είναι διπλά , δηλαδή σημείων που d > ε
				KDTREEProsthesimidiplwn(vrn,vtr, e);
			}
			f = 1;
		}
		cout << "Plithos vertices: " << vtr.size() << endl;
		getPointCloud2(vtr, pcM);
		auto t5 = std::chrono::high_resolution_clock::now();
		auto dt3 = std::chrono::duration_cast<std::chrono::seconds>(t5 - t4).count();
		cout << "XRONOS IPOLOGISMOY OLIKOY NEFOYS M: " << dt3 << " seconds" << endl;
		q1 = 0;
	}




	///////////    ERWTIMA 6

	
	if ((er6 == true) && (q1 == 1)) {
		int M = 0;
		int f = 0;
		float e = 0;
		vrn.clear();
		vsobel.clear();
		float dt = 0;
		for (int i = Ni; i < Nj; i++) {
			/// Convert the image to Grey
			cv::Mat grey;
			cv::cvtColor(I[i], grey, CV_BGR2GRAY);
			cv::Mat sobelx;
			cv::Sobel(grey, sobelx, CV_8U, 1, 0, 3);
			double minVal, maxVal;
			minMaxLoc(sobelx, &minVal, &maxVal); //find minimum and maximum intensities
			/// Sobel Image
			cv::Mat sobelxc;
			sobelx.convertTo(sobelxc, CV_8U, 255.0 / (maxVal - minVal), -minVal * 255.0 / (maxVal - minVal));
			/// Threshold
			cv::Mat threshold;
			threshold = sobelxc > 10;
			/// Display Sobel and Threshold Images
			cv::imshow("Sobel", sobelxc);
			cv::imshow("Threshold", threshold);

			/// Image i  -> pc1
			PCSobel(F[i], vertices1, threshold);
			getPointCloud2(vertices1, pc1);
			/// Image i+1 -> pc2
			PCSobel(F[i + 1], vertices2, threshold);
			getPointCloud2(vertices2, pc2);
			auto t6 = std::chrono::high_resolution_clock::now();
			float edistance = 0;
			std::vector<vec> P;
			if (f == 0) {     //   Na ekteleitai mono tin prwti fora - First Pass
				for (auto& v : vertices1)
				{
					vsobel.push_back(v);
				}
				/// Euresi plisiesterwn vertices tou PC Ni ws pros to PC Ni+1
				PlisiesterosKDTREE(vertices2, vertices1, P);
				/// Eksagwgi PC plisiesterwn vertices
				vvr::PointCloud PNearest;
				getPointCloud2(P, PNearest);
				/// Rotation - Translation - Μετατόπιση και Περιστροφή καθενός Pointcloud(i+1) προς τα πλησιέστερα που ανήκουν στο Pointcloud(i)
				RotTran RT;
				RT = rotran(pc2, PNearest);   // Μετατόπιση και Περιστροφή καθενός Pointcloud ως προς
				if (pc1.size() >= pc2.size()) {			 //   με τα πλησιέστερα του προηγούμενου
					M = pc2.size();
				}
				if (pc1.size() < pc2.size()) {
					M = pc1.size();
				}					                   
				for (int j = 0; j < M; j++) {
					float vv1 = (RT.R(0, 0)*pc2[j].x + RT.R(0, 1)*pc2[j].y + RT.R(0, 2)*pc2[j].z + RT.t(0));
					float vv2 = (RT.R(1, 0)*pc2[j].x + RT.R(1, 1)*pc2[j].y + RT.R(1, 2)*pc2[j].z + RT.t(1));
					float vv3 = (RT.R(2, 0)*pc2[j].x + RT.R(2, 1)*pc2[j].y + RT.R(2, 2)*pc2[j].z + RT.t(2));
					Point3D pe(vv1, vv2, vv3);
					/// Αποστάσεις των πλησιέστερων από τα vertices του Ni+1 για υπολογισμό της ανοχής
					edistance = edistance + (vertices2[j].DistanceSq(P[j]));
					/// Εισαγωγή στο vrn των μετατοπισμένων σημείων
					vrn.push_back(vec(vv1, vv2, vv3));   //  Metatopismena vertices
				}
				/// Ως ανοχή e θεωρείται η μέση απόσταση των πλησιέστερων vertices tou PC Ni ws pros to PC Ni+1 από τα vertices του Ni+1
				/// Προσοχή!! Υπολογίζεται μόνο στο first Pass και στις επόμενες επαναλήψεις παραμένει το ίδιο 
				e = edistance / float(M);
				cout << "Anoxi ε: " << e << endl;
				/// Εισαγωγή στο vtr των ευθυγραμμισμένων σημείων του Νi+1 ως προς το Νi που δεν είναι διπλά , δηλαδή σημείων που d > ε
				KDTREEProsthesimidiplwn(vrn, vsobel, e);
			}
			else {
				/// Euresi plisiesterwn vertices tou euthigrammismenou PC tis proigoumenis epanalipsis i ws pros to PC i+1
				PlisiesterosKDTREE(vertices2, vrn, P);
				/// Eksagwgi PC plisiesterwn vertices
				vvr::PointCloud PNearest;
				getPointCloud2(P, PNearest);
				/// Rotation - Translation - Μετατόπιση και Περιστροφή καθενός Pointcloud(i+1) προς τα πλησιέστερα που ανήκουν στο Pointcloud(i)
				RotTran RT;
				RT = rotran(pc2, PNearest);   // Μετατόπιση και Περιστροφή καθενός Pointcloud ως προς
				int M = pc2.size();                      //   με τα πλησιέστερα του προηγούμενου ευθυγραμμισμένου
				vrn.clear();
				for (int j = 0; j < M; j++) {
					float vv1 = (RT.R(0, 0)*pc2[j].x + RT.R(0, 1)*pc2[j].y + RT.R(0, 2)*pc2[j].z + RT.t(0));
					float vv2 = (RT.R(1, 0)*pc2[j].x + RT.R(1, 1)*pc2[j].y + RT.R(1, 2)*pc2[j].z + RT.t(1));
					float vv3 = (RT.R(2, 0)*pc2[j].x + RT.R(2, 1)*pc2[j].y + RT.R(2, 2)*pc2[j].z + RT.t(2));
					Point3D pe(vv1, vv2, vv3);
					/// Εισαγωγή στο vrn των μετατοπισμένων σημείων
					vrn.push_back(vec(vv1, vv2, vv3));
				}
				/// Εισαγωγή στο vtr των ευθυγραμμισμένων σημείων του i+1 ως προς το i που δεν είναι διπλά , δηλαδή σημείων που d > ε
				KDTREEProsthesimidiplwn(vrn, vsobel, e);
			}
			auto t7 = std::chrono::high_resolution_clock::now();
			auto dt4 = std::chrono::duration_cast<std::chrono::seconds>(t7 - t6).count();
			dt = dt + dt4;
			f = 1;
		}
		cout << "Plithos vertices: " << vsobel.size() << endl;
		getPointCloud2(vsobel, pcS);
		cout << "BELTIWMENOS XRONOS IPOLOGISMOY OLIKOY NEFOYS M ME SOBEL: " << dt << " seconds" << endl;
		q1 = 0;
	}

	



	//////////           DRAW           ///////////

	
	if (er1 == true) {         // ERWTIMA 1
		drawPointCloud(pointclouds[N], vvr::Colour::green);
	}
	

	if (er2 == true) {         // ERWTIMA 2
		drawPointCloud(pointclouds[N], vvr::Colour::red);
		if (N != 0) {
			drawPointCloud(pointclouds[N - 1], vvr::Colour::yellow);
		}
		for (auto& lt : ls1)
		{
			lt.draw();
		}
	}
	
	if (er4 == true) {         // ERWTIMA 4
		//drawPointCloud(pnew1, vvr::Colour::black);
		drawPointCloud(pcM, vvr::Colour::magenta);
	}

	
	///   Erwtima Aii me KDTREES

	if ((er5 == true) && (er8 == true)) {       // ERWTIMA 5 A
		for (auto& lt : ls)
		{
			lt.draw();
		}
		drawPointCloud(pointclouds[N5], vvr::Colour::black);
		if (N5 != 0) {
			drawPointCloud(pointclouds[N5 - 1], vvr::Colour::magenta);
		}
	}

	///   Erwtimata 3,4 me KDTREES

	if ((er5 == true) && (er9 == true)) {         // ERWTIMA 5 B
		drawPointCloud(pcM, vvr::Colour::cyan);
	}	

		
	if (er6 == true) {         // ERWTIMA 6
		getPointCloud2(vsobel, pcS);
		drawPointCloud(pcS, vvr::Colour::darkGreen);
		drawPointCloud(pcsobel1, vvr::Colour::grey);
	}

	if (er7 == true) {         // ERWTIMA 7
		for (int i = 0; i < tris.size(); i++) {
			vvr::Triangle &t = tris[i];
			Triangle3D td(
				t.v1().x, t.v1().y, t.v1().z,
				t.v2().x, t.v2().y, t.v2().z,
				t.v3().x, t.v3().y, t.v3().z,
				Colour::yellow);
			td.draw();
		}

	}

	// Display Images and update OpenCV
	// IMPORTANT: OpenCV finds open windows by name, e.g. "Color Image"
	
	cv::imshow("Color Image", m_color_image);
	cv::imshow("Depth Image", m_depth_image_norm);
	cv::imshow("Depthsssss Image", m_depth);
	
	
	// Following is necessary to update CV's loop
	cv::waitKey(10);

}



int main(int argc, char* argv[])
{
	
    try {
        return vvr::mainLoop(argc, argv, new PointCloudScene);
    }
    catch (std::string exc) {
        cerr << exc << endl;
        return 1;
    }
    catch (...)
    {
        cerr << "Unknown exception" << endl;
        return 1;
    }
}

void PointCloudScene::getPointCloud(const vvr::Mesh& src, vvr::PointCloud& dst, const vvr::Colour& color)
{
	std::vector<vec> vertices = src.getVertices();
	dst.clear();

	for(auto& pt : vertices)
	{
		dst.push_back(Point3D(pt.x, pt.y, pt.z, color));
	}
}

void PointCloudScene::drawPointCloud(const vvr::PointCloud& cloud, const vvr::Colour& color)
{
	
	for (auto& pt : cloud) {
		Point3D ptn(pt.x, pt.y, pt.z, color);
		ptn.draw();
	}
}

void PC(const cv::Mat& depth, std::vector<vec> &vertices)
{

	float fovx = 62;
	float fovy = 48.6;
	float	focal_length_constant_x = 1/(2 * tan((fovx*(pi) / 180) / 2)); 
	float	focal_length_constant_y = 1/(2 * tan((fovy*(pi) / 180) / 2));

	float ResX = depth.cols;
	float ResY = depth.rows;
	float f_x = ResX * focal_length_constant_x;
	float f_y = ResY * focal_length_constant_y;
	float Vo = depth.cols / 2;
	float Uo = depth.rows / 2;
	for (int i = 0; i < depth.cols; i = i + D) {           // Δειγματοληψία i = i + D
		for (int j = 0; j < depth.rows; j = j + D) {          // Δειγματοληψία j = j + D
			float d = depth.at<unsigned char>(j, i);
			float Z = d;
			float X = (i - Vo) * (Z / f_x);
			float Y = (j - Uo) * (Z / f_y);
			//vertices[j*depth.cols + i] = (vec(X, -Y, -Z));
			vertices.push_back(vec(X, -Y, -Z));
		}
	}
}

void PointCloudScene::sliderChanged(int slider_id, float v)
{
	switch (slider_id)
	{
	case 0:
		if (v != 0) {
			l0 = N;
			N = (v)*(178 + 1);
			echo(N);
			if (N != l0) {
				q0 = 1;
			}
		}
		break;
	case 1:
		if (v != 0) {
			l1 = Ni;
			Ni = (v)*(178 + 1);
			echo(Ni);
			if (Ni != l1) {
				q1 = 1;
			}
		}
		break;
	case 2:
		if (v != 0) {
			l2 = Nj;
			Nj = (v)*(178 + 1);
			echo(Nj);
			if (Nj != l2) {
				q1 = 1;
			}
		}
		break;
	case 3:
		if (v != 0) {
			l3 = N5;
			N5 = (v)*(178 + 1);
			echo(N5);
			if (N5 != l3) {
				q0ii = 1;
			}
		}
		break;
	}
}


void PointCloudScene::getPointCloud2(std::vector<vec> &vertices, vvr::PointCloud& dst)
{
	dst.clear();

	for (auto& pt : vertices)
	{
		dst.push_back(Point3D(pt.x, pt.y, pt.z));
	}
}

void NearestNeighboorKDtree(const vec& test_pt, const KDNode* root, const KDNode* &nn, float& best_dist)
{
	if (!root) return;

	int axis = root->axis;
	double distance = test_pt.ptr()[axis] - root->split_point.ptr()[axis];
	double best_dist_t = test_pt.DistanceSq(root->split_point);

	if (best_dist_t < best_dist) {
		nn = root;
		best_dist = best_dist_t;
	}
	bool isRight = distance > 0;
	NearestNeighboorKDtree(test_pt, isRight ? root->child_right : root->child_left, nn, best_dist);

	if (distance*distance < best_dist) {
		NearestNeighboorKDtree(test_pt, isRight ? root->child_left : root->child_right, nn, best_dist);
	}
}

void PCSobel(const cv::Mat& depth, std::vector<vec> &vertices, const cv::Mat& thresh)
{

	float fovx = 62;
	float fovy = 48.6;
	float	focal_length_constant_x = 1/(2 * tan((fovx*(pi) / 180) / 2));  
	float	focal_length_constant_y = 1/(2 * tan((fovy*(pi) / 180) / 2));

	float ResX = depth.cols;
	float ResY = depth.rows;
	float f_x = ResX * focal_length_constant_x;
	float f_y = ResY * focal_length_constant_y;
	float Vo = depth.cols / 2;
	float Uo = depth.rows / 2;
	for (int i = 0; i < depth.cols; i = i + D)   {      // Δειγματοληψία i = i + D
		for (int j = 0; j < depth.rows; j = j + D) {          // Δειγματοληψία j = j + D
			if (thresh.at<unsigned char>(j, i) == 255) {    // Pairnoume to epipedo 255 , diladi ta aspra (Ta 0 einai ta maura)
				float d = depth.at<unsigned char>(j, i);
				float Z = d;
				float X = (i - Vo) * (Z / f_x);
				float Y = (j - Uo) * (Z / f_y);
				//vertices[j*depth.cols + i] = (vec(X, -Y, -Z));
				vertices.push_back(vec(X, -Y, -Z));
			}
		}
	}
}

void Plisiesteros(std::vector<vec> vertices1, std::vector<vec> vertices2, std::vector<vec> &P, float SinApostasi, float MesiApostasi) {

	int N1 = vertices1.size();
	int N2 = vertices2.size();
	vec p;
	for (int i = 0; i < N1; i++) {      // N1
		float minDistance = std::numeric_limits<float>::max();
		int k = 0;
		for (int j = 0; j < N2; j++) {              //    N2
			float pointDistance = vertices1[i].DistanceSq(vertices2[j]);
			if (pointDistance < minDistance) {
				minDistance = pointDistance;
				p.x = vertices2[j].x;
				p.y = vertices2[j].y;
				p.z = vertices2[j].z;
			}		
		}
		P.push_back(p);
		SinApostasi = SinApostasi + minDistance; 
	}
	MesiApostasi = SinApostasi / N1;

}

void PointCloudScene::PlisiesterosKDTREE(std::vector<vec> vertices1, std::vector<vec> vertices2, std::vector<vec> &P) {
	m_KDTree = new KDTree(vertices2);
	int S = vertices1.size();
	float best_dist;
	vec nn;
	for (int i = 0; i < S; i++) {
		best_dist = std::numeric_limits<float>::max();
		const KDNode *nearest = NULL;
		NearestNeighboorKDtree(vertices1[i], m_KDTree->root(), nearest, best_dist);
		if (nearest) nn = nearest->split_point;
		P.push_back(nn);
	}
	delete m_KDTree;
}

void PointCloudScene::KDTREEProsthesimidiplwn(std::vector<vec> vertices1, std::vector<vec> &vertices2, float e) {
	m_KDTree = new KDTree(vertices2);
	int S = vertices1.size();
	float best_dist;
	vec nn;
	for (int i = 0; i < S; i++) {
		best_dist = std::numeric_limits<float>::max();
		const KDNode *nearest = NULL;
		NearestNeighboorKDtree(vertices1[i], m_KDTree->root(), nearest, best_dist);
		if (nearest) nn = nearest->split_point;
		float Distance = nn.DistanceSq(vertices1[i]);
		if (Distance > e) {
			vertices2.push_back(vertices1[i]);
		}
	}
	delete m_KDTree;
}