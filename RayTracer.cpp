
#include "RayTracer.h"
#include "GlobalData.h"
#include <fstream>

#include <omp.h>


/*	Plan of action
* 1. Shoot ray of light with given wavelength
*		Sampled from some distribution
* 2. Trace the ray through the scene
*		Interacting with glass use Fresnel's equation for intensity of each color
*		Wavelength bends depending on the index of refraction of that wavelength for this medium
*			- Using the dispersion function of hydrated silica (if i can find it)
* 3. (Theoretically) At each step of intersections, find nearby photons, compare their directions
*		Collect photons that have similar directions
*		With these compute interference, new intesities (possibly combine photons/rays)
*
* 
* 
* We shoot a ray, a ray is some coherent bundle of photons of the same frequency
* Each photon has its own frequency and phase while the group of photons defines the 'amplitude' or 'intensity' of the ray
* 
* 
* 
* Intensity at a point is proportiona to the square of the average amplitude of the wave
*/

/* 
* Other notes look at the txt file
*/


#define SRGB_ALPHA 0.055

inline double linear_to_srgb(double x) {
	if (x <= 0.0031308)
		return 12.92 * x;
	else
		return (1 + SRGB_ALPHA) * (pow(x, 1 / 2.4) - SRGB_ALPHA);
}

inline double pw_gauss(double x, double mu, double s1, double s2) {
	return (x < mu) ? exp(-0.5 * pow(x - mu, 2) / (s1 * s1)) : exp(-0.5 * pow(x - mu, 2) / (s2 * s2));
}

// // Operates on wavelengths in NANOMETERS (visible 380-780)
inline Vec3 WavelengthToRGB(double lambda) {
	double x, y, z, r, g, b;
	x = 1.056 * pw_gauss(lambda, 599.8, 37.9, 31.0) + 0.362 * pw_gauss(lambda, 442.0, 16.0, 26.7) -
		0.065 * pw_gauss(lambda, 501.1, 20.4, 26.2);
	y = 0.821 * pw_gauss(lambda, 568.8, 46.9, 40.5) + 0.286 * pw_gauss(lambda, 530.9, 16.3, 31.1);
	z = 1.217 * pw_gauss(lambda, 437.0, 11.8, 36.0) + 0.681 * pw_gauss(lambda, 459.0, 26.0, 13.8);

	r = 3.2406 * x - 1.5372 * y - 0.4986 * z;
	g = -0.9689 * x + 1.8758 * y + 0.0415 * z;
	b = 0.0557 * x - 0.2040 * y + 1.0570 * z;

	r = linear_to_srgb(r);
	g = linear_to_srgb(g);
	b = linear_to_srgb(b);

	r = fmin(fmax(r, 0.0), 1.0);
	g = fmin(fmax(g, 0.0), 1.0);
	b = fmin(fmax(b, 0.0), 1.0);

	return Vec3(r, g, b);
}

// Generates a Triclinic packed volume of spheres
void RayTracer::GenerateTriclinicPack(const Vec3& offset, const int width, const int height, const int layers) {
	double rad = 0.350, gappedRad = rad + 0.01;
	double xstep = 2.0 * gappedRad, ystep = sqrt(8.0 / 3.0) * gappedRad, zstep = sqrt(3.0) * gappedRad;
	double dx = gappedRad, dz = sqrt(1.0 / 3.0) * gappedRad;

	Vec3 sc = Vec3(0.1, 0.9, 0.4);

	for (int k = 0; k < layers; k++) {			// y-axis
		for (int j = 0; j < height; j++) {		// z-axis
			for (int i = 0; i < width; i++) {	// x-axis
				shapes.push_back(new Sphere(Vec3(offset.x + i * xstep + j * gappedRad + k * dx, 
					offset.y + k * ystep, 
					offset.z + j * zstep + k * dz), rad, sc, Vec3(0)));
			}
		}
	}
}


double sceneTheta = PI / 2.0;
void RayTracer::createScene() {
	rayTree.clear();
	while (!shapes.empty()) {
		Primitive* temp = shapes.back();
		shapes.pop_back();
		delete temp;
	}

	while (!lights.empty()) {
		Light* temp = lights.back();
		lights.pop_back();
		delete temp;
	}

	shapes.push_back(new Sphere(Vec3(-5, 0, -5), 0.1, Vec3(1.00, 0.1, 0.1)));


	// Opal Plane
	double ps = 2.0;	// Plane size
	double dy = 0;		// Vertical shift

	for (int i = 0; i < GLOBAL_data->numOpalLayers; i++) {
		Plane* opal = new Plane(Vec3(-ps, dy, -ps), Vec3(-ps, dy, ps), Vec3(ps, dy, ps), Vec3(ps, dy, -ps), Vec3(0.1, 0.1, 0.9));
		opal->SetMatType(Opal);
		shapes.push_back(opal);

		// Have the layers be physically split by micrometers
		// Shift by the width of the silica spheres
		dy -= GLOBAL_data->opalSphereSize / 1000.0;

		//Plane* water = new Plane(Vec3(-ps, dy, -ps), Vec3(-ps, dy, ps), Vec3(ps, dy, ps), Vec3(ps, dy, -ps), Vec3(0.5, 0.8, 0.2));
		//water->SetMatType(Water);
		//shapes.push_back(water);

		//// Shift by the width of the water
		//dy -= GLOBAL_data->waterRatio * (GLOBAL_data->opalSphereSize / 1000.0);
	}
	// Bottom reflective layer
	Plane* refl = new Plane(Vec3(-ps, dy, -ps), Vec3(-ps, dy, ps), Vec3(ps, dy, ps), Vec3(ps, dy, -ps), Vec3(0.2, 0.2, 0.2));
	refl->SetMatType(Reflector);
	shapes.push_back(refl);


	// Creating the light and observing plane
	double theta = sceneTheta, R = 4, R2 = 3;
	double sideLen = 0.5, sideLen2 = 0.05;
	Vec3 radial1 = Vec3(R * cos(theta), R * sin(theta), 0);
	Vec3 z = Vec3(0, 0, 1);
	Vec3 inPlane = radial1.cross(z);
	inPlane.normalize();

	z *= sideLen;
	inPlane *= sideLen;
	// Light Plane
	Plane t1 = Plane(radial1 + inPlane + z, radial1 - inPlane + z, radial1 - inPlane - z, radial1 + inPlane - z, Vec3(0.2, 0.2, 0.2));

	Vec3 radial2 = Vec3(R2 * cos(theta), R2 * sin(theta), 0);
	z = Vec3(0, 0, 1);
	inPlane = radial2.cross(z);
	inPlane.normalize();

	z *= sideLen2;
	inPlane *= sideLen2;
	// Light Cone Plane
	Vec3 c1 = radial2 + inPlane + z,
		c2 = radial2 - inPlane + z,
		c3 = radial2 - inPlane - z,
		c4 = radial2 + inPlane - z;
	Plane t2 = Plane(c1, c2, c3, c4, Vec3(0.2, 0.2, 0.2));
	lights.push_back(new Light(t1, t2));


	// Observing plane
	theta = PI - theta;
	sideLen = 4.0 * sideLen;
	radial1.x *= -1;
	z = Vec3(0, 0, 1);
	inPlane = radial1.cross(z);
	inPlane.normalize();

	z *= sideLen;
	inPlane *= sideLen;
	Plane * coll = new Plane(radial1 + inPlane + z, radial1 - inPlane + z, radial1 - inPlane - z, radial1 + inPlane - z, Vec3(0.2, 0.2, 0.2));
	coll->SetMatType(Collector);
	shapes.push_back(coll);
	collector = coll;
}

RayTracer::RayTracer() {
    texture = 0;

	dx = 0;
	dy = 0;
	dz = 0;
	moveScale = 0.001;

	camPos = Vec3(0.0, 0.0, 0.0);
	camDir = Vec3(0.0, 0.0, -1.0);
	camUp = Vec3(0.0, 1.0, 0.0);


	// position, radius, surface color, emission color, transparency, reflectivity
	//spheres.push_back(new Sphere(Vec3(0, -10004, -20), 10000, Vec3(0.2), 0.0, 0.0));

	//shapes.push_back(new Sphere(Vec3(3, 0, -15), 2, Vec3(1.00, 0.1, 0.1), Vec3(0), 0.95, 0.65));
	//shapes.push_back(new Sphere(Vec3(1, -1, -18), 1, Vec3(1.0, 1.0, 1.0), Vec3(0), 0.9, 0.9));
	//shapes.push_back(new Sphere(Vec3(-2, 2, -15), 2, Vec3(0.1, 0.1, 1.0), Vec3(0), 0.5, 0.05));
	//shapes.push_back(new Sphere(Vec3(-4, 3, -18), 1, Vec3(0.1, 1.0, 0.1), Vec3(0), 0.7, 0.3));
	//shapes.push_back(new Sphere(Vec3(-4, 0, -25), 1, Vec3(1.00, 0.1, 0.1), Vec3(0), 0.95, 0.65));
	//shapes.push_back(new Sphere(Vec3(-1, 1, -25), 2, Vec3(1.0, 1.0, 1.0), Vec3(0), 0.0, 0.0));
	//shapes.push_back(new Sphere(Vec3(2, 2, -25), 1, Vec3(0.1, 0.1, 1.0), Vec3(0), 0.5, 0.05));
	//shapes.push_back(new Sphere(Vec3(5, 3, -25), 2, Vec3(0.1, 1.0, 0.1), Vec3(0), 0.7, 0.3));
	////shapes.push_back(new Sphere(Vec3(-10, 20, 0), 3, Vec3(0), Vec3(3), 0, 0));
	////shapes.push_back(new Sphere(Vec3(0, 10, 0), 3, Vec3(0), Vec3(1), 0, 0));

	//shapes.push_back(new Plane(Vec3(1, -4, 1), Vec3(1, -4, 7), Vec3(4, -4, 4), Vec3(4, -4, 1), Vec3(0.5, 0.1, 0.5), Vec3(0), 0, 0));
	
	//shapes.push_back(new Plane(Vec3(0, 0, 0), Vec3(2, 0, 0), Vec3(2, 2, 0), Vec3(0, 2, 0), Vec3(0.1, 0.1, 1.0), Vec3(0), 0, 0));
	//shapes.push_back(new Plane(Vec3(0, 0, -1), Vec3(0, 2, -1), Vec3(2, 2, -1), Vec3(2, 0, -1), Vec3(0.5, 0.1, 0.5), Vec3(0), 0, 0));

	//shapes.push_back(new Plane(Vec3(-5, -5, -2), Vec3(-5, 5, -2), Vec3(5, 5, -2), Vec3(5, -5, -2), Vec3(0.5, 0.1, 0.5), Vec3(0), 0, 0));

	//GenerateTriclinicPack(Vec3(0, 0, -4), 4, 2, 4);


	//shapes.push_back(new Sphere(Vec3(0, 0, 0), 0.1, Vec3(0.9, 0.1, 0.1)));
}

RayTracer::~RayTracer() {
	while (!shapes.empty()) {
		Primitive* temp = shapes.back();
		shapes.pop_back();
		delete temp;
	}

	while (!lights.empty()) {
		Light* temp = lights.back();
		lights.pop_back();
		delete temp;
	}
}

/*
	Shooting rays and intersecting and drawing stuff :)
*/


/*
	DRAWING STUFF
*/
Ray RayTracer::ComputeReflectance(const Vec3& dir, const Vec3& normal, const double& n1, const double& n2, double& reflectance, bool method) {
	double cosi = dir.dot(-normal), n12 = n1 / n2;
	double sint = n12 * sqrt(fmax(0.0, 1.0 - cosi * cosi));

	// Check for total internal reflection
	if (sint >= 1.0) {
		reflectance = 1.0;
		return Ray(dir + (2.0 * cosi) * normal, Vec3());
	}

	double cost = sqrt(fmax(0.0, 1.0 - sint * sint));
	if (method == 0) {
		// Schlick's approximation
		double R0 = pow((n1 - n2) / (n1 + n2), 2);
		//return R0 + (1.0 - R0) * pow(1.0 - cosi, 5);
		reflectance = R0 + (1.0 - R0) * pow(1.0 - cosi, 5);
	}
	else {
		// Fresnel averaging method
		double n1cosi = n1 * cosi, n1cost = n1 * cost,
			n2cosi = n2 * cosi, n2cost = n2 * cost;
		double Rs = (n1cosi - n2cost) / (n1cosi + n2cost);
		double Rp = (n1cost - n2cosi) / (n1cost + n2cosi);
		//return (Rs * Rs + Rp * Rp) / 2.0;
		reflectance = (Rs * Rs + Rp * Rp) / 2.0;
	}
	// Ray used to store the reflection vector and the transmitted vector
	return Ray((dir + (2.0 * cosi) * normal).normalize(), (n12 * dir + (n12 * cosi - cost) * normal).normalize());
}

// Operates on wavelengths in MICROMETERS 10^-6 (visible 0.380-0.780)
double RayTracer::SilicaDispersionFunction(double lambda) {
	double l2 = lambda * lambda;
	double n21 = (A1 * l2 / (l2 - A2)) + (B1 * l2 / (l2 - B2)) + (C1 * l2 / (l2 - C2));
	return sqrt(n21 + 1);
}

// Operates on wavelengths in MICROMETERS 10^-6 (visible 0.380-0.780)
double RayTracer::WaterDispersionFunction(double lambda) {
	double l2 = lambda * lambda;
	double n21 = (X1 * l2 / (l2 - X2)) + (Y1 * l2 / (l2 - Y2));
	return sqrt(n21 + 1);
}


// We go into this method once the ray has hit an 'Opal' material
// Since the opal is the a singlular plane we 
void RayTracer::OpalTraceRay(std::vector<Ray>& contactPhotons, const Vec3& orig, Ray* ray, Hit& hit, int bounces, int layer, bool visible) {

	bool intersects = false;
	// Check intersection with each shape
	for (Primitive* p : shapes) {
		if (p->Intersect(ray, hit)) intersects = true;
	}

	// TECHNICALLY if it intersects and the surface is reflective
	Vec3 hitPoint = ray->PointAtParameter(hit.GetT());

	if (visible) {
		Vec3 color = ray->energy * WavelengthToRGB(ray->wavelength);
		rayTree.push_back(Line(ray->origin, hitPoint, color));
	}
	if (intersects && bounces > 0 && hit.GetMaterial()->GetMatType() != Collector) {
		double n1 = 1, refl;
		double n2;
		if (hit.GetMaterial()->GetMatType() == Water)
			n2 = WaterDispersionFunction(ray->wavelength / 1000.0);
		else
			n2 = SilicaDispersionFunction(ray->wavelength / 1000.0);
		if (!hit.GetFrontFace()) {
			double t = n1;
			n1 = n2;
			n2 = t;
		}


		double lambda = ray->wavelength / 1000.0;
		// Calculate lambda inside the medium
		lambda = lambda / n2;
		double modulo = fmod(hit.GetT(), lambda);
		ray->AddPhase(modulo / lambda);


		// Temp stores the reflected direction and transmission direction
		Ray temp = ComputeReflectance(ray->dir, hit.GetNormal(), n1, n2, refl, 0);
		double reflEnergy = ray->energy * refl;
		double refrEnergy = ray->energy * (1.0 - refl);

		Hit reflectedHit;

		int nextLayer = layer;
		if (hit.GetFrontFace()) nextLayer--;
		else nextLayer++;
		if (reflEnergy > GLOBAL_data->energyThresh)
			OpalTraceRay(contactPhotons, orig, new Ray(hitPoint + EPSILON * temp.origin, temp.origin, ray->wavelength, reflEnergy, ray->phase), reflectedHit, bounces - 1, nextLayer, visible);
		
		if (refl < 1.0 && hit.GetMaterial()->GetMatType() != Reflector && refrEnergy > GLOBAL_data->energyThresh) {
			Hit refractedHit;

			nextLayer = layer;
			if (hit.GetFrontFace()) nextLayer++;
			else nextLayer--;

			if (hit.GetMaterial()->GetMatType() == Opal && layer == 0) {
				//Vec3 sp = (hitPoint + EPSILON * temp.dir) - orig;
				//sp.x *= pow(10, -6);
				//sp.z *= pow(10, -6);
				//TraceRay(contactPhotons, new Ray(orig, temp.dir, ray->wavelength, refrEnergy, ray->phase), refractedHit, bounces - 1, visible);
				TraceRay(contactPhotons, new Ray(hitPoint + EPSILON * temp.dir, temp.dir, ray->wavelength, refrEnergy, ray->phase), refractedHit, bounces - 1, visible);
			}
			else
				OpalTraceRay(contactPhotons, orig, new Ray(hitPoint + EPSILON * temp.dir, temp.dir, ray->wavelength, refrEnergy, ray->phase), refractedHit, bounces - 1, nextLayer, visible);
		}
	}

	delete ray;
}


// Recursive function traces each ray
void RayTracer::TraceRay(std::vector<Ray>& contactPhotons, Ray* ray, Hit& hit, int bounces, bool visible) {

	bool intersects = false;
	// Check intersection with each shape
	for (Primitive* p : shapes) {
		if (p->Intersect(ray, hit)) intersects = true;
	}

	// TECHNICALLY if it intersects and the surface is reflective
	Vec3 hitPoint = ray->PointAtParameter(hit.GetT());

	if (visible) {
		Vec3 color = ray->energy * WavelengthToRGB(ray->wavelength);
		rayTree.push_back(Line(ray->origin, hitPoint, color));
	}
	if (intersects && bounces > 0) {
		//double n1 = 0.5, refl;
		//double n2 = SilicaDispersionFunction(ray->wavelength / 1000.0);
		////n2 = SilicaDispersionFunction(ray->wavelength / 1000.0);
		////double n2 = (ray->wavelength - 0.380) / 0.400 + 1;
		//if (!hit.GetFrontFace()) {
		//	double t = n1;
		//	n1 = n2;
		//	n2 = t;
		//}
		//if (hit.GetMaterial()->GetMatType() != Collector) {
		//	Ray temp = ComputeReflectance(ray->dir, hit.GetNormal(), n1, n2, refl, 0);
		//	Hit reflectedHit;
		//	TraceRay(new Ray(hitPoint + EPSILON * temp.origin, temp.origin, ray->wavelength, ray->amplitude, ray->phase), reflectedHit, bounces - 1, visible);
		//	if (refl < 1.0) {
		//		Hit refractedHit;
		//		TraceRay(new Ray(hitPoint + EPSILON * temp.dir, temp.dir, ray->wavelength, ray->amplitude, ray->phase), refractedHit, bounces - 1, visible);
		//	}
		//}

		// If we've hit something then add to the phase shift
		double lambda = ray->wavelength * pow(10, -9);
		double modulo = fmod(hit.GetT(), lambda);
		ray->AddPhase(modulo / lambda);

		if (hit.GetMaterial()->GetMatType() != Collector) {
			double n1 = 1.0, refl;
			double n2 = SilicaDispersionFunction(ray->wavelength / 1000.0);
			if (!hit.GetFrontFace()) {
				double t = n1;
				n1 = n2;
				n2 = t;
			}
			// Temp stores the reflected direction and transmission direction
			Ray temp = ComputeReflectance(ray->dir, hit.GetNormal(), n1, n2, refl, 0);
			double reflEnergy = ray->energy * refl;
			double refrEnergy = ray->energy * (1.0 - refl);

			Hit reflectedHit;
			if (reflEnergy > GLOBAL_data->energyThresh)
				TraceRay(contactPhotons, new Ray(hitPoint + EPSILON * temp.origin, temp.origin, ray->wavelength, reflEnergy, ray->phase), reflectedHit, bounces - 1, visible);
			if (refl < 1.0 && refrEnergy > GLOBAL_data->energyThresh) {
				Hit refractedHit;
				OpalTraceRay(contactPhotons, hitPoint + EPSILON * temp.origin, new Ray(hitPoint + EPSILON * temp.dir, temp.dir, ray->wavelength, refrEnergy, ray->phase), refractedHit, bounces - 1, 1, visible);
			}
		}

		if (hit.GetMaterial()->GetMatType() == Collector) {
			contactPhotons.push_back(Ray(ray));
		}
	}

	delete ray;
	return;
}

// Traces rays from the current viewing direction
void RayTracer::TraceViewRay() {
	rayTree.clear();
	Hit hit;

	/*for (double lm = 380; lm <= 780; lm += 10) {
		TraceRay(Ray(camPos, camDir, lm, 1.0), hit, GLOBAL_data->numBounces, true);
	}*/
	photons.clear();
	//TraceRay(new Ray(camPos, camDir, 650, 1.0, 0), hit, GLOBAL_data->numBounces, true);

	std::vector<Ray> contactPhotons;
	TraceRay(contactPhotons, new Ray(camPos, camDir, 650, 1.0, 0), hit, GLOBAL_data->numBounces, true);
	std::cout << "Total photons: " << photons.size() << std::endl;
}
// Traces rays from the provided light source
void RayTracer::TraceLightRays() {
	rayTree.clear();
	photons.clear();

	//std::ofstream file;
	//file.open("C:\\MyStuff\\RPI_STUFF\\Project_OPAL\\PythonSim\\test2.txt");

	//for (; GLOBAL_data->numOpalLayers <= 10; GLOBAL_data->numOpalLayers += 1) {
	//for (; sceneTheta >= PI/16.0; sceneTheta -= PI/32.0) {
		createScene();

		for (Light* l : lights) {
			Vec3 center = l->getCenter();
			Vec3 normal = l->getNormal();

			Vec3 dx = l->proj.p2 - l->proj.p1;
			Vec3 dy = l->proj.p4 - l->proj.p1;
			Vec3 randp;

			for (double lm = 380; lm <= 780; lm += 1) {
				std::vector<Ray> contactPhotons;

				randp = l->proj.p1 + GLOBAL_data->rand() * dx + GLOBAL_data->rand() * dy;
				Hit hit;
				TraceRay(contactPhotons, new Ray(center + normal * 0.001, (randp - center).normalize(), lm, 1.0, 0), hit, GLOBAL_data->numBounces, true);

				double totalEnergy = 0;
				for (Ray phot : contactPhotons) {
					totalEnergy += phot.energy * sin(phot.phase * 2.0 * PI);
				}
				totalEnergy = fabs(totalEnergy);
				//file << GLOBAL_data->numOpalLayers << " " << lm << " " << totalEnergy << "\n";

				if (fmod(lm, 50) < 1)
					std::cout << GLOBAL_data->numOpalLayers << " On lambda = " << lm << std::endl;
			}
		}
	//}
	//file.close();
}
void RayTracer::DrawRayTree() {
	for (Line& r : rayTree) {
		r.Draw();
	}
}

void RayTracer::RayTrace() {
	glColor3f(1.0f, 0, 0);
	glBegin(GL_LINES);

	Vec3 perp = camDir.cross(camUp);
	perp.normalize();

	rayTree.clear();
	double angle = tan(PI * 0.5 * GLOBAL_data->fov / 180.0);
	for (int y = 0; y < GLOBAL_data->height; y += 50) {
		for (int x = 0; x < GLOBAL_data->width; x += 50) {
			double xx = (2 * ((x + 0.5) * GLOBAL_data->invW) - 1) * angle * GLOBAL_data->aspectRatio;
			double yy = (1 - 2 * ((y + 0.5) * GLOBAL_data->invH)) * angle;

			Vec3 dir(xx, yy, -1);
			dir.normalize();
			//dir *= 5.0;

			dir = dir.rotateAboutV(camUp, totYaw);
			dir = dir.rotateAboutV(perp, totPitch);

			Hit hit;
			//TraceRay(Ray(Vec3(0, 0, 0), dir), hit, 1, true);
			glVertex3f(0.0f, 0.0f, 0.0f);
			glVertex3f(dir.x, dir.y, dir.z);

		}
	}
	glEnd();
}


Vec3 randVec() {
	Vec3 tmp;
	while (true) {
		tmp = Vec3(2 * GLOBAL_data->rand() - 1,  // random real in [-1,1]
			2 * GLOBAL_data->rand() - 1,  // random real in [-1,1]
			2 * GLOBAL_data->rand() - 1); // random real in [-1,1]
		if (tmp.norm() < 1) break;
	}
	tmp.normalize();
	return tmp;
}

// Function to set the camera movement based on keyboard input
void RayTracer::SetCamMovement(char key, bool down) {
	int modif = glutGetModifiers();

	if (key == 'e' && down) {
		using namespace std;
		//for (double i = 0.210; i <= 1.0; i += 0.010) {
		//	cout << i << " = " << DispersionFunction(i) << endl;
		//}


		rayTree.clear();
		int count = 1000000;
		Ray dir(Vec3(2.5, 2.0, -7), Vec3(0, 0, 1));
		// create display list

		dlist = glGenLists(1);
		glNewList(dlist, GL_COMPILE);
		glBegin(GL_POINTS);

		for (int i = 0; i < count; ++i) {
			Ray randRay(dir.origin, (dir.dir + randVec()).normalize());

			bool intersects = false;
			Hit hit;
			for (Primitive* p : shapes) {
				if (p->Intersect(&randRay, hit)) intersects = true;
			}

			if (intersects) {
				Vec3 hitPoint = randRay.PointAtParameter(hit.GetT());
				glColor4d(1.0, 0, 0, 1.0);
				glVertex4d(hitPoint.x, hitPoint.y, hitPoint.z, 1.0);
			}
		}
		glEnd();
		glEndList();
	}

	if (key == 't' && down) {
		TraceViewRay();
		//TraceLightRays();
		return;
	}

	if (key == 'p' && down) {
		sceneTheta -= PI / 256.0;
		createScene();
	}
	if (key == 'o' && down) {
		sceneTheta += PI / 256.0;
		createScene();
	}
	if (key == 'i' && down) {
		TraceLightRays();
	}
	if (key == 'b' && down) {
		createScene();
	}

	switch (key) {
	case 'w':		// forward
		dz = down ? 1.0 : 0.0;
		break;
	case 's':		// backward
		dz = down ? -1.0 : 0.0;
		break;
	case 'd':		// strafe right
		dx = down ? 1.0 : 0.0;
		break;
	case 'a':		// strafe left
		dx = down ? -1.0 : 0.0;
		break;
	case ' ':
		if (modif == GLUT_ACTIVE_SHIFT)
			dy = down ? -1.0 : 0.0;
		else
			dy = down ? 1.0 : 0.0;
		break;
	default:
		break;
	}
}

void RayTracer::RotateCam(int dyaw, int dpitch) {
	double pitch = dpitch / 200.0 * TO_RADS;
	double yaw = dyaw / 200.0 * TO_RADS;

	if (pitch != 0 || yaw != 0) {
		camDir = camDir.rotateAboutV(camUp, yaw);

		Vec3 perp = camDir.cross(camUp);
		perp.normalize();
		camDir = camDir.rotateAboutV(perp, pitch);

		totPitch += pitch;
		totYaw = fmod(totYaw + yaw, 2.0 * PI);
	}
}

// Helper lookAt function to use custom vector
void LookAt(Vec3 pos, Vec3 dir, Vec3 up) {
	gluLookAt(pos.x, pos.y, pos.z,
		pos.x + dir.x, pos.y + dir.y, pos.z + dir.z,
		up.x, up.y, up.z);
}
void RayTracer::Display() {
	// Clear and reset display settings
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Update the camera positions based on any user input
	Vec3 perp = camDir.cross(camUp);
	perp.normalize();
	camPos += (dz * moveScale * camDir);
	camPos += (dy * moveScale * camUp);
	camPos += (dx * moveScale * perp);

	// Update the view matrix so that we are looking from the new camera position
	glLoadIdentity();
	LookAt(camPos, camDir, camUp);


	// Draw all our spheres
	for (Primitive* s : shapes) {
		s->Draw();
	}
	for (Light* s : lights) {
		s->Draw();
	}


	//RayTrace();
	DrawRayTree();

	//glRotatef(GLOBAL_data->angle, 0.0f, 0.0f, 1.0f);
	//glRotatef(GLOBAL_data->angle * 0.5f, 1.0f, 0.0f, 0.0f);
	//double lambda = (GLOBAL_data->angle / (64 * PI)) * 400.0 + 380.0;
	//Vec3 v = WavelengthToRGB(lambda);
	//glColor3ub((int)(v.x * 255), (int)(v.y * 255), (v.z * 255));
	//glBegin(GL_TRIANGLES);
	//glVertex3f(0, -2, 3.0);
	//glVertex3f(4, 0.0, 3.0);
	//glVertex3f(2, 2, 3.0);
	//glEnd();
	//if (GLOBAL_data->angle > 64 * PI)
	//	GLOBAL_data->angle = 0;
	//else
	//	GLOBAL_data->angle += 0.001f;


	glCallList(dlist);

	//glFlush();
	glutSwapBuffers();
}