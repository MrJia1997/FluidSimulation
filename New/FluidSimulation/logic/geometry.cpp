#include "constants.h"
#include "geometry.h"

#include <cmath>
#include <algorithm>

CubeContainer::CubeContainer():
    particleMap(100, triple_hash),
    xmin(X_MIN), xmax(X_MAX),
    ymin(Y_MIN), ymax(Y_MAX),
    zmin(Z_MIN), zmax(Z_MAX)
{
    bucketSize = KERNEL_H;
}

CubeContainer::~CubeContainer()
{
    // do nothing
}

void CubeContainer::clear()
{
	particleMap.clear();
}

IntTriple CubeContainer::indexOfPosition(const Particle& p) {
	IntTriple index;
	index.x = (int)std::floor(p.predPosition[0] / bucketSize);
	index.y = (int)std::floor(p.predPosition[1] / bucketSize);
	index.z = (int)std::floor(p.predPosition[2] / bucketSize);
	return index;
}

IntTriple CubeContainer::indexOfPosition(const double& x,const double& y, const double& z) {
    IntTriple index;
    index.x = (int)std::floor(x / bucketSize);
    index.y = (int)std::floor(y / bucketSize);
    index.z = (int)std::floor(z / bucketSize);
    return index;
}

void CubeContainer::add(int i, const Particle& p)
{
	IntTriple index = indexOfPosition(p);
	
    if (particleMap.count(index) > 0)
    {
		//Then the bucket already exists, so just add this to the bin
		particleMap[index].push_back(i);
	}
	// Otherwise we need to create the bucket
    else
    {
        particleMap[index] = std::vector<int>();
		particleMap[index].push_back(i);
	}
}

bool closeEnough(const Particle& p1, const Particle& p2)
{
    QVector3D diff = p1.predPosition - p2.predPosition;
    if (std::fabs(diff[0]) >= KERNEL_H || std::fabs(diff[1]) >= KERNEL_H || std::fabs(diff[2]) >= KERNEL_H)
    {
		return false;
	}
    return (diff.length() <= KERNEL_H);
}

void CubeContainer::findNeighbors(Particle& p_i, std::vector<Particle>& particles)
{
	IntTriple index = indexOfPosition(p_i);
	p_i.neighbors.clear();

	// Look in all neighboring grid cells
    for (int dx = -1; dx <= 1; dx++)
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            for (int dz = -1; dz <= 1; dz++)
            {
				IntTriple neighborIndex = index.addOffset(dx, dy, dz);
				if (particleMap.count(neighborIndex) == 0) continue;
                for (auto& j: particleMap[neighborIndex])
                {
                    if (closeEnough(p_i, particles[j]))
                    {
						p_i.neighbors.push_back(j);
					}
				}
			}
		}
	}
}

void CubeContainer::handleCollision(Particle& particle)
{
    // double lossFactor = 0.6;

    QVector3D afterDelta = particle.predPosition;

    /*if (afterDelta.x() < xmin)
    {
        particle.predPosition[0] = 2 * xmin - particle.predPosition.x();
        //particle.velocity.x() = (-1.0) * lossFactor * particle.velocity.x();
    }
    else if (afterDelta.x() > xmax)
    {
        particle.predPosition[0] = 2 * xmax - particle.predPosition.x();
        //particle.velocity.x() = (-1.0) * lossFactor * particle.velocity.x();
    }

    if (afterDelta.y() < ymin)
    {
        particle.predPosition[1] = 2 * ymin - particle.predPosition.y();
        //particle.velocity.y() = (-1.0) * lossFactor * particle.velocity.y();
    }
    else if (afterDelta.y() > ymax)
    {
        particle.predPosition[1] = 2 * ymax - particle.predPosition.y();
        //particle.velocity.y() = (-1.0) * lossFactor * particle.velocity.y();
    }

    if (afterDelta.z() < zmin)
    {
        particle.predPosition[2] = 2 * zmin - particle.predPosition.z();
        //particle.velocity.z() = (-1.0) * lossFactor * particle.velocity.z();
    }
    else if (afterDelta.z() > zmax)
    {
        particle.predPosition[2] = 2 * zmax - particle.predPosition.z();
        //particle.velocity.z() = (-1.0) * lossFactor * particle.velocity.z();
    }*/

    /*if (afterDelta.x() < xmin)
    {
        particle.predPosition[0] = xmin;
        //particle.velocity.x() = (-1.0) * lossFactor * particle.velocity.x();
    }
    else if (afterDelta.x() > xmax)
    {
        particle.predPosition[0] = xmax;
        //particle.velocity.x() = (-1.0) * lossFactor * particle.velocity.x();
    }

    if (afterDelta.y() < ymin)
    {
        particle.predPosition[1] = ymin;
        //particle.velocity.y() = (-1.0) * lossFactor * particle.velocity.y();
    }
    else if (afterDelta.y() > ymax)
    {
        particle.predPosition[1] = ymax;
        //particle.velocity.y() = (-1.0) * lossFactor * particle.velocity.y();
    }

    if (afterDelta.z() < zmin)
    {
        particle.predPosition[2] = zmin;
        //particle.velocity.z() = (-1.0) * lossFactor * particle.velocity.z();
    }
    else if (afterDelta.z() > zmax)
    {
        particle.predPosition[2] = zmax;
        //particle.velocity.z() = (-1.0) * lossFactor * particle.velocity.z();
    }*/

    double diss = 4.0;

    if (afterDelta.x() < xmin)
    {
        particle.predPosition[0] = xmin + (xmin - particle.predPosition.x()) / diss;
        //particle.velocity.x() = (-1.0) * lossFactor * particle.velocity.x();
    }
    else if (afterDelta.x() > xmax)
    {
        particle.predPosition[0] = xmax + (xmax - particle.predPosition.x()) / diss;
        //particle.velocity.x() = (-1.0) * lossFactor * particle.velocity.x();
    }

    if (afterDelta.y() < ymin)
    {
        particle.predPosition[1] = ymin + (ymin - particle.predPosition.y()) / diss;
        //particle.velocity.y() = (-1.0) * lossFactor * particle.velocity.y();
    }
    else if (afterDelta.y() > ymax)
    {
        particle.predPosition[1] = ymax + (ymax - particle.predPosition.y()) / diss;
        //particle.velocity.y() = (-1.0) * lossFactor * particle.velocity.y();
    }

    if (afterDelta.z() < zmin)
    {
        particle.predPosition[2] = zmin + (zmin - particle.predPosition.z()) / diss;
        //particle.velocity.z() = (-1.0) * lossFactor * particle.velocity.z();
    }
    else if (afterDelta.z() > zmax)
    {
        particle.predPosition[2] = zmax + (zmax - particle.predPosition.z()) / diss;
        //particle.velocity.z() = (-1.0) * lossFactor * particle.velocity.z();
    }
}
