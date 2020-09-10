#include <VVRScene/canvas.h>
#include <VVRScene/mesh.h>
#include <MathGeoLib.h>
#include <opencv2/highgui.hpp>
#include <Eigen\Dense>



namespace vvr {
	typedef std::vector<vvr::Point3D> PointCloud;
	std::vector<LineSeg3D> ls1;
	std::vector<LineSeg3D> ls2;
	std::vector<LineSeg3D> ls;
}

typedef std::vector<vec> VecVector;

struct KDNode
{
	vec split_point;     // to simeio tou root
	int axis;
	int level;
	AABB aabb;
	KDNode *child_left;
	KDNode *child_right;
	KDNode() : child_left(NULL), child_right(NULL) {}
	~KDNode() { delete child_left; delete child_right; }
};

class KDTree
{
public:
	KDTree(VecVector &pts);
	~KDTree();
	std::vector<KDNode*> getNodesOfLevel(int level);
	int depth() const { return m_depth; }
	const KDNode* root() const { return m_root; }
	const VecVector &pts;

private:
	static int makeNode(KDNode *node, VecVector &pts, const int level);
	static void getNodesOfLevel(KDNode *node, std::vector<KDNode*> &nodes, int level);

private:
	KDNode *m_root;
	int m_depth;
};

struct VecComparator {
	unsigned axis;
	VecComparator(unsigned axis) : axis(axis % 3) {}
	virtual inline bool operator() (const vec& v1, const vec& v2) {
		return (v1.ptr()[axis] < v2.ptr()[axis]);
	}
};


class PointCloudScene : public vvr::Scene
{
	// Simple bit-controlled enumerator
	enum
	{
		FLAG_SHOW_SOLID = 1 << 0,
		FLAG_SHOW_NORMALS = 1 << 1,
		FLAG_SHOW_WIRE = 1 << 2,
		FLAG_SHOW_AXES = 1 << 3,
		FLAG_SHOW_POINTS = 1 << 4,
		FLAG_SHOW_PLANE = 1 << 5,
		FLAG_NUM
	};
public:
	PointCloudScene();
	const char* getName() const { return "PointCloud Scene"; }
	void keyEvent(unsigned char key, bool up, int modif) override;
	void arrowEvent(vvr::ArrowDir dir, int modif) override;
	void PointCloudScene::sliderChanged(int slider_id, float v);


private:
	void draw() override;
	void reset() override;
	void resize() override;

	// Helper method to get and draw PointCloud from Mesh
	void getPointCloud(const vvr::Mesh& src, vvr::PointCloud& dst, const vvr::Colour& color);
	void drawPointCloud(const vvr::PointCloud& cloud, const vvr::Colour& color);
	void getPointCloud2(std::vector<vec> &vertices, vvr::PointCloud& dst);
	void PlisiesterosKDTREE(std::vector<vec> vertices1, std::vector<vec> vertices2, std::vector<vec> &P);
	void PointCloudScene::KDTREEProsthesimidiplwn(std::vector<vec> vertices1, std::vector<vec> &vertices2, float e);


	// Members
	int m_style_flag;
	float m_plane_d;
	vvr::Canvas2D m_canvas;
	vvr::Colour m_point_color;
	vvr::Mesh m_model_original, m_model;
	math::Plane m_plane;
	vvr::PointCloud m_model_points;
	vvr::PointCloud m_model_points1;
	vvr::PointCloud m_model_points2;
	vvr::PointCloud pc1;
	vvr::PointCloud pc2;
	std::vector<cv::Mat> F;
	std::vector<cv::Mat> I;
	std::vector<cv::Mat> FI;
	std::vector<vvr::PointCloud> pointclouds;
	std::vector<vvr::PointCloud> pointcloudnew;
	std::vector<vec> vertices1;
	std::vector<vec> vertices2;
	std::vector<vec> ver;
	std::vector<vec> ver7;
	std::vector<vec> v1;
	std::vector<vec> v2;
	std::vector<vec> vtr;
	std::vector<vec> vrn;
	std::vector<vec> vsobel;
	std::vector<vec> vplisiesteros;
	std::vector<vec> modelVerts;
	std::vector<std::vector <vec>> vall;
	vvr::PointCloud pnew1;
	vvr::PointCloud pcM;
	vvr::PointCloud pcS;
	vvr::PointCloud pcsobel1;
	vvr::PointCloud pcsobel2;
	KDTree *m_KDTree;
	std::vector<vvr::Triangle> tris;


	// example depth image
	cv::Mat m_depth_image, m_color_image, m_depth_image_norm, m_depth_image_new, m_depth_image_norm_previous, m_depth;

};
void PC(const cv::Mat& depth, std::vector<vec> &vertices);
void NearestNeighboorKDtree(const vec& test_pt, const KDNode* root, const KDNode* &nn, float& best_dist);
void Plisiesteros(std::vector<vec> vertices1, std::vector<vec> vertices2, std::vector<vec> &P, float SinApostasi, float MesiApostasi);
void PCSobel(const cv::Mat& depth, std::vector<vec> &vertices, const cv::Mat& thresh);