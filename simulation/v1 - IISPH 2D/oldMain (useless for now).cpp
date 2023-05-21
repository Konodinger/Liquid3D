// ----------------------------------------------------------------------------
// main.cpp
//
//  Created on: Fri Jan 22 20:45:07 2021
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: SPH simulator (DO NOT DISTRIBUTE!)
//
// Copyright 2021-2023 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
// ----------------------------------------------------------------------------

#define _USE_MATH_DEFINES

#include <stdio.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <math.h>
//#include <omp.h>

#include "Vector.hpp"
#include "IISPH_solver.cpp"

using namespace std;


// window parameters
/*GLFWwindow *gWindow = nullptr;
int gWindowWidth = 1024;
int gWindowHeight = 768;*/

// timer
float timer = 0.f;
bool gSolverStop = false;

// global options
/*bool gPause = true;
bool gSaveFile = false;
bool gSaveMovie = false;
bool gShowGrid = true;
bool gShowVel = false;
int gSavedCnt = 0;

const int kViewScale = 20;*/

SphSolver gSolver(0.08, 0.5, 1e3, Vec2f(0, -9.8), 0.01, 7.0, 0.5, 0.5, 1.);

void printHelp()
{
  cout <<
    "> Help:" << endl <<
    "    Keyboard commands:" << endl <<
    "    * H: print this help" << endl <<
    "    * P: toggle simulation" << endl <<
    "    * G: toggle grid rendering" << endl <<
    "    * V: toggle velocity rendering" << endl <<
    "    * S: save current frame into a file" << endl <<
    "    * B: start or end registering frames into a video" << endl <<
    "    * Q: quit the program" << endl;
}

// Executed each time the window is resized. Adjust the aspect ratio and the rendering viewport to the current window.
/*void windowSizeCallback(GLFWwindow *window, int width, int height)
{
  ///
}*/

// Executed each time a key is entered.
/*void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
  if(action == GLFW_PRESS && key == GLFW_KEY_H) {
    printHelp();
  } else if(action == GLFW_PRESS && key == GLFW_KEY_S) {
    gSaveFile = !gSaveFile;
  } else if(action == GLFW_PRESS && key == GLFW_KEY_B) {
    if (gSaveMovie) {
      cout << "Recording video..." << endl;
      const char* cmd = "ffmpeg -y -r 60 -f image2 -i capture/s%04d.tga -vcodec libx264 -crf 25  -pix_fmt yuv420p output.mp4";
      ffmpeg = _popen(cmd, "wb");
      _pclose(ffmpeg);
      ffmpeg = _popen("rm -r -f capture", "wb");
    _pclose(ffmpeg);
    } else {
      ffmpeg = _popen("md capture", "wb");
      _pclose(ffmpeg);
      gSavedCnt = 0;
    }
    gSaveMovie = !gSaveMovie;
  } else if(action == GLFW_PRESS && key == GLFW_KEY_G) {
    gShowGrid = !gShowGrid;
  } else if(action == GLFW_PRESS && key == GLFW_KEY_V) {
    gShowVel = !gShowVel;
  } else if(action == GLFW_PRESS && key == GLFW_KEY_P) {
    gSolverStop = !gSolverStop;
    if(!gSolverStop)
      gAppTimerLastClockTime = static_cast<float>(glfwGetTime());
  } else if(action == GLFW_PRESS && key == GLFW_KEY_Q) {
    glfwSetWindowShouldClose(window, true);
  }
}*/

void init()
{
  gSolver.initScene(48, 32, 16, 16);
}

// The main rendering call
void render()
{
  ///

  // grid guides
  ///

  // render particles
  ///

  // velocity
  ///

  /*if (gSaveMovie) {
    stringstream fpath;
    fpath << "capture/s" << setw(4) << setfill('0') << gSavedCnt++ << ".tga";

    const short int w = gWindowWidth;
    const short int h = gWindowHeight;
    vector<int> buf(w*h*3, 0);
    glReadPixels(0, 0, w, h, GL_BGR, GL_UNSIGNED_BYTE, &(buf[0]));

    FILE *out = fopen(fpath.str().c_str(), "wb");
    short TGAhead[] = {0, 2, 0, 0, 0, 0, w, h, 24};
    fwrite(&TGAhead, sizeof(TGAhead), 1, out);
    fwrite(&(buf[0]), 3*w*h, 1, out);
    fclose(out);
  }*/

  /*if(gSaveFile) {
    stringstream fpath;
    fpath << "s" << setw(4) << setfill('0') << gSavedCnt++ << ".tga";

    cout << "Saving file " << fpath.str() << " ... " << flush;
    const short int w = gWindowWidth;
    const short int h = gWindowHeight;
    vector<int> buf(w*h*3, 0);
    glReadPixels(0, 0, w, h, GL_BGR, GL_UNSIGNED_BYTE, &(buf[0]));

    FILE *out = fopen(fpath.str().c_str(), "wb");
    short TGAhead[] = {0, 2, 0, 0, 0, 0, w, h, 24};
    fwrite(&TGAhead, sizeof(TGAhead), 1, out);
    fwrite(&(buf[0]), 3*w*h, 1, out);
    fclose(out);
    gSaveFile = false;

    cout << "Done" << endl;
  }*/
}

int main(int argc, char **argv)
{

  //omp_set_num_threads( 12 );

  init();
  while(true) {
    if (!gSolverStop){
      for(int i=0; i<10; ++i) gSolver.update();
      timer += 0.1f;
    }
    render();
  }
  cout << " > Quit" << endl;
  return EXIT_SUCCESS;
}
