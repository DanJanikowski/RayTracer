
#include "GlobalData.h"
#include "RayTracer.h"


GlobalData* GLOBAL_data;


GlobalData::GlobalData(int argc, char** argv) {
	angle = 0.0f;

	raytracer = new RayTracer();

	GLOBAL_data = this;
}

GlobalData::~GlobalData() {
	delete raytracer;
}