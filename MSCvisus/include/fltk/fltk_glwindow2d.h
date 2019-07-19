#ifndef FLTK_GLWINDOW2D_H
#define FLTK_GLWINDOW2D_H

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>

class fltkGlWindow2D : public Fl_Gl_Window {
public:

	float mScaleFactor;
	float mTranslateX;
	float mTranslateY;
	bool mFirstRender;

	bool mInvertBackground;
	bool mDrawBoundingBox;

	float extents[2];

	void FixViewport(int W,int H) {
		glLoadIdentity();
		glViewport(0,0,W,H);
		glOrtho(-W,W,-H,H,-1,1);
	}

	virtual void draw_rest() {}


	void draw() {
		if (mFirstRender) { 
			FixViewport(w(), h()); 
			mFirstRender=false; 
		}      

		// Clear screen to bg color
		if (mInvertBackground) {
			glClearColor(0,0,0, 0.0);
		} else {
			glClearColor(1,1,1, 1.0);
		}
		glClear(GL_COLOR_BUFFER_BIT);

		// set up scene transformations
		glPushMatrix();
		glScalef(mScaleFactor, mScaleFactor, 1.0);
		glTranslatef(2.0*mTranslateX,2.0*mTranslateY, 0.0);

		if (mDrawBoundingBox) {
			if (mInvertBackground) {
				glColor3f(1,1,1);
			} else {
				glColor3f(0,0,0);
			}
			glBegin(GL_LINE_LOOP);
			glVertex2f(0,0);
			glVertex2f(extents[0], 0);
			glVertex2f(extents[0], extents[1]);
			glVertex2f(0, extents[1]);
		}

		// THIS IS OVERLOADED FUNCTION
		draw_rest();

		glPopMatrix();
	}
	// HANDLE WINDOW RESIZING
	void resize(int X,int Y,int W,int H) {
		Fl_Gl_Window::resize(X,Y,W,H);
		FixViewport(W,H);
		redraw();
	}

	int oldx;
	int oldy;
	int pushedbutton;

	int handle(int e) {
		int ret = Fl_Gl_Window::handle(e);

		switch(e) {

		case FL_PUSH:
			oldx = Fl::event_x();
			oldy = Fl::event_y();
			pushedbutton = Fl::event_button();
			return 1;
		case FL_DRAG:
			int newx = Fl::event_x();
			int newy = Fl::event_y();

			if (pushedbutton == 3) {
				float diffy = 0.05 * (oldy - newy);
				mScaleFactor += diffy;
				if (mScaleFactor < 0.1) mScaleFactor = 0.01;
			} else if (pushedbutton == 1) {
				float diffx = (newx - oldx)/mScaleFactor;
				float diffy = (oldy - newy)/mScaleFactor;
				mTranslateX += diffx;
				mTranslateY += diffy;
			}

			oldx = newx;
			oldy = newy;
			redraw();
		}

		return ret;
	}

public:

	// OPENGL WINDOW CONSTRUCTOR
	fltkGlWindow2D(int extentX, int extentY, int X,int Y,int W,int H,const char*L=0) : Fl_Gl_Window(X,Y,W,H,L) {
		extents[0] = extentX;
		extents[1] = extentY;
		mFirstRender = true;
		mInvertBackground = true;
		mDrawBoundingBox = true;
		mScaleFactor = 1.0; 
		mTranslateX = 0.0;
		mTranslateY = 0.0;
		mIdArray = NULL;
		select_o_matic = false;
		mSelectedTexturer = 0;
		end();
	}

};



#endif