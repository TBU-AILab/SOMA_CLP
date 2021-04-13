#ifndef KMEANS_H
#define KMEANS_H

#include <vector>
#include <cmath>
#include <algorithm>

class Point {

private:
	int pointId, clusterId;
	int dimensions;
	std::vector<double> values;
	double CF;

public:
	Point(int id, double* in, double CFin, int dim) {
		dimensions = dim;
		pointId = id;
		CF = CFin;
		for (int i = 0; i < dimensions; i++) {
			values.push_back(in[i]);
		}
		clusterId = 0; //Initially not assigned to any cluster
	}

	int getDimensions() {
		return dimensions;
	}

	int getCluster() {
		return clusterId;
	}

	int getID() {
		return pointId;
	}

	void setCluster(int val) {
		clusterId = val;
	}

	double getVal(int pos) {
		return values[pos];
	}

	double getCF() {
		return CF;
	}

	void setCF(double CFin) {
		CF = CFin;
	}

};

class Cluster {

private:
	int clusterId;
	std::vector<double> centroid;
	std::vector<Point> points;

public:
	Cluster(int clusterId, Point centroid) {
		this->clusterId = clusterId;
		for (int i = 0; i < centroid.getDimensions(); i++) {
			this->centroid.push_back(centroid.getVal(i));
		}
		this->addPoint(centroid);
	}

	void addPoint(Point p);

	bool removePoint(int pointId);

	int getId();

	Point getPoint(int pos);

	int getSize();

	double getCentroidByPos(int pos);

	void setCentroidByPos(int pos, double val);
};

class KMeans {
private:
	int K, iters, dimensions, total_points;
	std::vector<Cluster> clusters;

	int getNearestClusterId(Point point);

public:
	KMeans(int K, int iterations) {
		this->K = K;
		this->iters = iterations;
	}

	void run(std::vector<Point>& all_points);
};

#endif // KMEANS_H