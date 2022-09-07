#include <iostream>
#include <functional>

#include "GlobalData.h"
#include "OpenGLCanvas.h"


int main(int argc, char** argv) {

	GlobalData data(argc, argv);

	OpenGLCanvas canvas(argc, argv);

	return 1;
}