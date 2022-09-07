
#include "OpenGLCanvas.h"
#include "RayTracer.h"


void display();
void reshape(int w, int h);
void keyboardPress(unsigned char key, int x, int y);
void keyboardRelease(unsigned char key, int x, int y);
void cameraMove(int x, int y);
void timer(int);

OpenGLCanvas::OpenGLCanvas(int& argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(GLOBAL_data->width, GLOBAL_data->height);
	glutInitWindowPosition(200, 200);
	glutCreateWindow("ProjectOPAL");

	glutSetCursor(GLUT_CURSOR_NONE);
	glutDisplayFunc(display);
	glutIdleFunc(display);
	glutReshapeFunc(reshape);
	glutTimerFunc(0, timer, 0);

	glutKeyboardFunc(keyboardPress);
	glutKeyboardUpFunc(keyboardRelease);
	glutPassiveMotionFunc(cameraMove);

	// OpenGL init
	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_CULL_FACE);


	////GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	////GLfloat mat_shininess[] = { 50.0 };
	//GLfloat light_position[] = { -5.0, 10.0, 0.0, 0.0 };
	//glClearColor(0.0, 0.0, 0.0, 0.0);
	//glShadeModel(GL_SMOOTH);
	////glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	////glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	//glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);


	timer(0);
	glutMainLoop();
}


void display() {
	GLOBAL_data->raytracer->Display();
}

void timer(int) {
	glutPostRedisplay();
	glutWarpPointer(GLOBAL_data->width / 2, GLOBAL_data->height / 2);
	glutTimerFunc(1000 / GLOBAL_data->FPS, timer, 0);
}

void reshape(int w, int h) {
	if (h == 0)
		h = 1;
	float ratio = 1.0 * w / float(h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, w, h);
	gluPerspective(GLOBAL_data->fov, ratio, GLOBAL_data->zNear, GLOBAL_data->zFar);
	glMatrixMode(GL_MODELVIEW);
}

void keyboardPress(unsigned char key, int x, int y) {
	GLOBAL_data->raytracer->SetCamMovement(key, true);

	switch (key) {
	case 27:
		exit(0);
		break;
	default:
		break;
	}
}

void keyboardRelease(unsigned char key, int x, int y) {
	GLOBAL_data->raytracer->SetCamMovement(key, false);
}

void cameraMove(int x, int y) {
	int dx, dy;
	dx = (GLOBAL_data->width / 2) - x;
	dy = (GLOBAL_data->height / 2) - y;

	GLOBAL_data->raytracer->RotateCam(dx, dy);
}