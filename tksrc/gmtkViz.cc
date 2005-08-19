// -*- C++ -*-

// Dots and dashes have been disabled because of the instability they
// induce. The wxWidgets documentation says that they may only be
// drawn with width one. I believe the problem is that when you print
// (and maybe scale, etc.) that width changes implicitly. There was a
// (perhaps) similar problem with fonts previously. The didn't scale
// the way they should have. That was solved by making a temporary
// font object based on the original font object's settings. See
// StructPage::draw(). Perhaps something similar can be done the
// pens. Instead of simply changing to the appropriate pen, create a
// new pen with the properties of that pen and switch to the new pen.

// the only parts of gmtk that gmtkViz needs to concern itself with
#include "GMTK_FileParser.h"
#include "GMTK_RVInfo.h"
// all of the wxWidgets headers for the things gmtkViz uses
#include <wx/wx.h>
#include <wx/cmdline.h>
#include <wx/colordlg.h>
#include <wx/dcbuffer.h>
#include <wx/ffile.h>
#include <wx/fontdlg.h>
#include <wx/gdicmn.h>
#include <wx/image.h>
#include <wx/list.h>
#include <wx/notebook.h>
#include <wx/print.h>
#include <wx/printdlg.h>
#include <wx/statline.h>
#include <wx/textdlg.h>
#include <wx/textfile.h>
#include <wx/tooltip.h>
#include <wx/tipwin.h>
#include <wx/textctrl.h>
#include <wx/spinctrl.h>
#include <wx/choice.h>
#include <wx/event.h>
#include <wx/font.h>
#include <cassert>
#include <cmath>
#include <vector>
#include <string.h>
#include <stack>
//needed for dirname and basename
#include <libgen.h>
//needed for stat
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <wx/filename.h>

//for .Xdefaults
#ifndef __WXMSW__
#include<X11/Xlib.h>
#include<X11/Xutil.h>
#endif

// Apprently these are needed in order to do anything with gmtk (even
// if you don't actually use them).
#include "GMTK_GMParms.h"
#include "rand.h"
#include "GMTK_ObservationMatrix.h"
GMParms GM_Parms;
RAND rnd(0);
ObservationMatrix globalObservationMatrix;

// Actually draw things 16x larger than we're going to display
// them. That way, if they want to zoom in 16x it still looks
// good. This should probably be a variable. Then we could just change
// it to whatever it needs to be to draw correctly.
#define ACTUAL_SCALE (1)
// The different zooming modes (in the menu) correspond to the
// following scales:
static const double gZoomMap[] = {
//	0.0625/ACTUAL_SCALE, // 1/16
//	0.125/ACTUAL_SCALE, // 1/8
//	0.25/ACTUAL_SCALE, // 1/4
	0.5/ACTUAL_SCALE, // 1/2
	0.70710678118654746/ACTUAL_SCALE, // 1/(2**(1/2))
	0.84089641525371461/ACTUAL_SCALE, // 1/(2**(1/4))
	0.91700404320467122/ACTUAL_SCALE, // 1/(2**(1/8))
	1.0/ACTUAL_SCALE, // 1
	1.0905077326652577/ACTUAL_SCALE, // 2**(1/8)
	1.189207115002721/ACTUAL_SCALE, // 2**(1/4)
	1.4142135623730951/ACTUAL_SCALE, // 2**(1/2)
	2.0/ACTUAL_SCALE, // 2
	4.0/ACTUAL_SCALE // 4
//	8.0/ACTUAL_SCALE, // 8
//	16.0/ACTUAL_SCALE // 16
};
//this is the index into gZoomMap which will give a Zoom of 1
#define ZOOM_1_INDEX 4

/// comparator for the maps
struct ltstr {
	bool operator()(const wxString s1, const wxString s2) const
	{
		return (s1 < s2);
	}
};

struct ltint {
	bool operator()(const int i1, const int i2) const
	{
		return (i1 < i2);
	}
};

//These store style data
//They are used so that we can save and restore styles of pens into a text file
//so the numeric style will be converted into text and when the text is read in
//it will be converted back to a number, this way if the numbers that represent
//the styles change within wxWidgets our files will still be valid, unless styles
//are dropped.  If new styles are added then just add them to the map where the
//others are added
//
//string -> int style (ie "wxSOLID" -> wxSOLID)
map<const wxString, int, ltstr> styleStringToIntMap;
//int style - > string (ie wxSOLID -> "wxSOLID")
map<const int, wxString, ltint> styleIntToStringMap;
//int style -> int index into pen_style_choices
map<const int, int, ltint> styleToStyleChoicesIndexMap;

//states of highlight
//off = not higlighed
//on = highlighted (for the children and parents)
//childern = highlighted and highlighting children
//parents = highlighted and highlighting parents
//both = highlighted and highlighting both parents and children
enum highlight_states {off,on,children,parents,both};


// some sizes of things
#define NODE_RADIUS (10*ACTUAL_SCALE)
#define GRID_SIZE (10*ACTUAL_SCALE)
#define NEW_CP_OFFSET (10*ACTUAL_SCALE)
#define ARROW_LEN (12*ACTUAL_SCALE)
#define ARROW_WID (4*ACTUAL_SCALE)

//forward declarations of things that GFrame needs
class wxgmtk1toManyDialog;
class	GmtkHelp;

// forward declarations of things StructPage needs
class NameTag;
class VizNode;
class VizArc;
class VizSep;
class Selectable;
class ControlPoint;

class gmtkPen;

static wxFont defaultFont;

//number of pens in a StructPage
#define numPens 14

//stores info about Pens
struct PenInfo {
	gmtkPen * pen;
	wxString nameabr;
	wxString namelong;
};

struct PenDefault_t {
	wxColour color;
	int style;
	int width;
};

map<const wxString, PenDefault_t *, ltstr> PenDefaultMap;

// this class subclasses wxPen and simply adds the ability to save "default" data
// and later see if the current pen setup matches that default.  This way we can 
// save only the non default pen data to a file
class gmtkPen : public wxPen {
	public: 
		gmtkPen(const wxColour& colour, int width, int style) : wxPen(colour,width,style){setDefaults();}
		gmtkPen(const wxString& colourName, int width, int style) : wxPen(colourName,width,style){setDefaults();}
//		gmtkPen(const wxBitmap& stipple, int width) : wxPen(wxBitmap,width){}
		gmtkPen(const wxPen& pen) : wxPen(pen){setDefaults();}
		gmtkPen(const gmtkPen& pen) : wxPen((wxPen)pen)
		{ default_style = pen.default_style; default_colour = pen.default_style; default_width = pen.default_width;}
		bool isDefault(){ return (default_colour == GetColour()) && (default_style == GetStyle()) && (default_width == GetWidth());}
		//save the current settings as the default
		void setDefaults(){ default_colour = GetColour(); default_style = GetStyle(); default_width = GetWidth();}
		//default setters and getters
		int GetDefaultStyle() {return default_style;}
		int GetDefaultWidth() {return default_width;}
		wxColour GetDefaultColour() {return default_colour;}
		void SetDefaultStyle(int new_style) {default_style = new_style;}
		void SetDefaultWidth(int new_width) {default_width = new_width;}
		void SetDefaultColour(wxColour new_colour) {default_colour = new_colour;}
	private:
		wxColour default_colour;
		int	default_style;
		int	default_width;
};


/*
 * This class creates a Dialog that allows users to customize all the pens 
 * at one time
 */

class GmtkPenSelectDialog : public wxDialog {
	public:
		GmtkPenSelectDialog( wxWindow *parent,
				wxWindowID id,
				const wxString &title,
				PenInfo * pen_info,
				int number_of_pens,
				const wxPoint& position = wxDefaultPosition,
				const wxSize& size = wxDefaultSize,
				long style = wxDEFAULT_DIALOG_STYLE );
	private:
		//event handlers
		void OnCancel( wxCommandEvent &event ){event.Skip();}
		void OnOk( wxCommandEvent &event ){event.Skip();}
		DECLARE_EVENT_TABLE()
};

/*
 * This class encapsulates the data and event handeler function
 * needed to react to the Color Customization dialog changes
 * in color
 */

class GmtkPenSelect_col_evthnd : public wxEvtHandler {
	public:
		GmtkPenSelect_col_evthnd(GmtkPenSelectDialog * parent_ptr, gmtkPen *
				pen_ptr, wxButton * colbtn_ptr) { pen = pen_ptr; colbutton =
			colbtn_ptr; parent = parent_ptr;}

		void OnEvent(wxCommandEvent& event)
		{
			GmtkPenSelect_col_evthnd * data =
				(GmtkPenSelect_col_evthnd*)event.m_callbackUserData;
			//pop up the get color dialog
			wxColour newColor = wxGetColourFromUser(data->parent,
					data->pen->GetColour());
			if (newColor.Ok()){
				//set the color of the pen and the button
				data->pen->SetColour(newColor);
				data->colbutton->SetBackgroundColour(newColor);
			}
			event.Skip();
		}
	private:
		gmtkPen * pen;
		wxButton * colbutton;
		GmtkPenSelectDialog * parent;
};

/*
 * This class encapsulates the data and event handeler function
 * needed to react to the Color Customization dialog changes
 * in width
 */

class GmtkPenSelect_width_evthnd : public wxEvtHandler {

	public:
		GmtkPenSelect_width_evthnd(GmtkPenSelectDialog * parent_ptr, gmtkPen *
				pen_ptr, wxSpinCtrl * wdth_sel_ptr) { pen = pen_ptr; width_select =
			wdth_sel_ptr; parent = parent_ptr; }

		void OnEvent(wxCommandEvent& event)
		{
			GmtkPenSelect_width_evthnd * data =
				(GmtkPenSelect_width_evthnd*)event.m_callbackUserData;
			int penWidth = data->width_select->GetValue();
#ifdef __WXMSW__
			//if we're using windows then we have to force the pen width to be 1
			//if the pen is non solid
			int penStyle = data->pen->GetStyle();
			if (penWidth != 1 && (penStyle == wxDOT || penStyle == wxLONG_DASH ||
					penStyle == wxSHORT_DASH || penStyle == wxDOT_DASH ||
					penStyle == wxUSER_DASH)) {
				wxLogMessage("Dots and dashes can only be used with pens of "
						"width 1, so you can't change the width until you change "
						"the style.");
				event.Skip();
				return;
			} 
#endif
			data->pen->SetWidth(penWidth);
			event.Skip();
		}
	private:
		gmtkPen * pen;
		wxSpinCtrl * width_select;
		GmtkPenSelectDialog * parent;
};

/*
 * This class encapsulates the data and event handeler function
 * needed to react to the Color Customization dialog changes
 * in style
 */
class GmtkPenSelect_style_evthnd : public wxEvtHandler {

	public:
		GmtkPenSelect_style_evthnd(GmtkPenSelectDialog * parent_ptr, gmtkPen *
				pen_ptr, wxChoice * style_choice_ptr) { pen = pen_ptr; style_choice
			= style_choice_ptr; parent = parent_ptr; }

		void OnEvent(wxCommandEvent& event)
		{
			//get the user data
			GmtkPenSelect_style_evthnd * data =
				(GmtkPenSelect_style_evthnd*)event.m_callbackUserData;
			int penStyle = data->style_choice->GetSelection();
			//set the pen style
			switch (penStyle) {
				case 0:
					data->pen->SetStyle(wxSOLID);
					break;
				case 1:
					data->pen->SetStyle(wxTRANSPARENT);
					break;
				case 2:
					data->pen->SetStyle(wxDOT);
					break;
				case 3:
					data->pen->SetStyle(wxLONG_DASH);
					break;
				case 4:
					data->pen->SetStyle(wxSHORT_DASH);
					break;
				case 5:
					data->pen->SetStyle(wxDOT_DASH);
					break;
					/*case 6:
					  data->pen->SetStyle(wxBDIAGONAL_HATCH);
					  break;
					  case 7:
					  data->pen->SetStyle(wxCROSSDIAG_HATCH);
					  break;
					  case 8:
					  data->pen->SetStyle(wxFDIAGONAL_HATCH);
					  break;
					  case 9:
					  data->pen->SetStyle(wxCROSS_HATCH);
					  break;
					  case 10:
					  data->pen->SetStyle(wxHORIZONTAL_HATCH);
					  break;
					  case 11:
					  data->pen->SetStyle(wxVERTICAL_HATCH);
					  break;*/
			}
#ifdef __WXMSW__
			//if we're using windows then we have to make sure that the
			//width is 1 if the pen is not solid or transparent
			if ( penStyle == wxDOT || penStyle == wxLONG_DASH ||
					penStyle == wxSHORT_DASH || penStyle == wxDOT_DASH ||
					penStyle == wxUSER_DASH ) {
				wxLogMessage("Dots and dashes can only be used with pens of "
						"width 1, so you can't change the width until you change "
						"the style.");
				data->pen->SetWidth(1);
			}
#endif
			event.Skip();
		}
	private:
		gmtkPen * pen;
		wxChoice * style_choice;
		GmtkPenSelectDialog * parent;
};


//the event table for the Pen Customization Dialog
BEGIN_EVENT_TABLE(GmtkPenSelectDialog,wxDialog)
    EVT_BUTTON( wxID_OK, GmtkPenSelectDialog::OnOk )
    EVT_BUTTON( wxID_CANCEL, GmtkPenSelectDialog::OnCancel )
END_EVENT_TABLE()

	//global variables that represent the number of pen styles
	//and the style names as the user sees them
	static const int n_pen_styles = 6;
	static const wxString pen_style_choices[] = { wxT("Solid"),
		wxT("Transparent"),
		wxT("Dotted"),
		wxT("Dashed (long dashes)"),
		wxT("Dashed (shart dashes)"),
		wxT("Dotted and Dashed") };

/*
 * GmtkPenSelectDialog Constructor, actually creates the dialog
 */

GmtkPenSelectDialog::GmtkPenSelectDialog( wxWindow *parent,
				wxWindowID id,
				const wxString &title,
				PenInfo * pen_info,
				int num_pens,
				const wxPoint& position,
				const wxSize& size,
				long style) : wxDialog(parent, id, title, position, size, style, "dialogBox")
{
	//allocate the color button data for event callbacks
	wxBoxSizer *topsizer = new wxBoxSizer( wxVERTICAL );
	wxFlexGridSizer * flex_sizer = new wxFlexGridSizer(4,5,5);
	//make labels
	flex_sizer->Add(new
			wxStaticText(this,-1,"Pen Name",wxDefaultPosition), 0,
			wxALL | wxALIGN_CENTER_VERTICAL | wxALIGN_CENTER_HORIZONTAL);
	flex_sizer->Add(new
			wxStaticText(this,-1,"Width",wxDefaultPosition), 0,
			wxALL | wxALIGN_CENTER_VERTICAL | wxALIGN_CENTER_HORIZONTAL);
	flex_sizer->Add(new
			wxStaticText(this,-1,"Style",wxDefaultPosition), 0,
			wxALL | wxALIGN_CENTER_VERTICAL | wxALIGN_CENTER_HORIZONTAL);
	flex_sizer->Add(new
			wxStaticText(this,-1,"Color",wxDefaultPosition), 0,
			wxALL | wxALIGN_CENTER_VERTICAL | wxALIGN_CENTER_HORIZONTAL);

	//add the pens
	for (int i = 0; i < num_pens; i++){
		//pen name
		flex_sizer->Add(new
				wxStaticText(this,-1,pen_info[i].namelong,wxDefaultPosition), 0,
				wxALL | wxALIGN_CENTER_VERTICAL);
		
		//width input
		wxSpinCtrl * width_select;
		flex_sizer->Add(width_select = new wxSpinCtrl(this,-1, wxEmptyString,
					wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS,  0, 16,
					pen_info[i].pen->GetWidth()), 0, wxALL |
				wxALIGN_CENTER_VERTICAL);
		//make an event handler object
		GmtkPenSelect_width_evthnd * width_evt_handler = new
			GmtkPenSelect_width_evthnd(this, pen_info[i].pen, width_select);
		//connect to the event handler
		width_select->Connect(-1, wxEVT_COMMAND_SPINCTRL_UPDATED,
				(wxObjectEventFunction) (wxEventFunction) (wxCommandEventFunction)
				&GmtkPenSelect_width_evthnd::OnEvent, (wxObject *)width_evt_handler);

		//style choice input
		wxChoice * style_choice;
		flex_sizer->Add(style_choice = new wxChoice(this, -1, wxDefaultPosition,
					wxDefaultSize, n_pen_styles, pen_style_choices), 0, wxALL |
				wxALIGN_CENTER_VERTICAL);
		style_choice->Select(styleToStyleChoicesIndexMap[pen_info[i].pen->GetStyle()]);
		//make an event handler object
		GmtkPenSelect_style_evthnd * style_evt_handler = new
			GmtkPenSelect_style_evthnd(this, pen_info[i].pen, style_choice);
		//connect to the event handler
		style_choice->Connect(-1, wxEVT_COMMAND_CHOICE_SELECTED,
				(wxObjectEventFunction) (wxEventFunction) (wxCommandEventFunction)
				&GmtkPenSelect_style_evthnd::OnEvent, (wxObject *)style_evt_handler);

		//color choice input
		wxButton * color_button;
		flex_sizer->Add(color_button = new wxButton(this, -1, wxT("Change"),
					wxDefaultPosition, wxDefaultSize), 0, wxALL |
				wxALIGN_CENTER_VERTICAL);
		color_button->SetBackgroundColour(pen_info[i].pen->GetColour());
		color_button->SetForegroundColour(*wxWHITE);
		//make an event handler object
		GmtkPenSelect_col_evthnd * col_evt_handler = new
			GmtkPenSelect_col_evthnd(this, pen_info[i].pen, color_button);
		//connect to the event handler
		color_button->Connect(-1, wxEVT_COMMAND_BUTTON_CLICKED,
				(wxObjectEventFunction) (wxEventFunction) (wxCommandEventFunction)
				&GmtkPenSelect_col_evthnd::OnEvent, (wxObject *)col_evt_handler);
	}

	//add the selections, align them along the center of the dialog box
	topsizer->Add(flex_sizer, 1, wxALIGN_CENTER | wxEXPAND | wxALL, 10);

	//create buttons and add them
	wxBoxSizer *button_sizer = new wxBoxSizer( wxHORIZONTAL );
	button_sizer->Add( new wxButton( this, wxID_OK, "Ok" ), 0, wxALL, 10 );
	button_sizer->Add( new wxButton( this, wxID_CANCEL, "Cancel" ), 0, wxALL, 10 );
	//add the buttons, align them along the center of the dialog box
	topsizer->Add( button_sizer, 0, wxALIGN_CENTER );

	SetAutoLayout( true );     // tell dialog to use sizer
	SetSizer( topsizer );      // actually set the sizer

	topsizer->Fit( this );            // set size to minimum size as calculated by the sizer
	topsizer->SetSizeHints( this );   // set size hints to honour mininum size
}

/// A tab with a scrolled area for displaying and manipulating a graph
/// associated with a structure file.
class StructPage: public wxScrolledWindow
{
	public:
		// constructor
		StructPage(wxWindow *parent, wxWindowID id,
				wxFrame *parentFrame, wxNotebook *parentNotebook,
				const wxString &file, bool old = true);
		// destructor
		virtual ~StructPage();

		//make sure that the data parsed ok
		bool parsedSucessfully() {return data_parsed;}

		// some event handlers
		void OnPaint( wxPaintEvent &event );
		void OnChar( wxKeyEvent &event );
		void OnMouseEvent( wxMouseEvent &event );
		void popUpNodeInfo(wxPoint);
		void OnCustomizePen();

		//functions called by event handlers
		void nextHighlightState(int node_index);
		void toggleNodeAndArcVisibility(int node_index);
		void selectAllInFrame(int frame);

		// pseudo event handlers: GFrame calls these.
		void Save( void );
		void SaveAs( void );
		bool RequestClose( void );
		// Set the status bar etc. when we're in the front
		void onComeForward( void );

		// Did everything parse successfully?
		bool Ready( void ) { return !gvpAborted; }

		// things related to selecting
		Selectable *itemAt( const wxPoint& pt );
		void setAllSelected( bool newSelected );
		void setEndpointsSelected( bool newSelected, RVInfo::rvParent );
		void setInOutArcsVisible( bool newVisible, RVInfo::rvParent rvId );
		void setInArcsHighlighed( bool newHighlighted, RVInfo::rvParent rvId );
		void toggleSelectedInRect( const wxRect& rect );
		void moveFrameSep( int i, int dx );
		void moveFrameNameTag( int i, int dx, int dy );
		void moveNodeNameTag( int i, int dx, int dy );
		void moveNode( int i, int dx, int dy );
		void moveControlPoint( int i, int j, int k, int dx, int dy );
		void moveSelected( int dx, int dy );
		void snapSelectedToGrid( void );

		void copyFrameLayout( int from, int to );
		void copyFrameLayout( void );
		void copyPartitionLayout( void );
		void copyArcLayout( int iFrom, int jFrom,
				int iTo, int jTo,
				bool backward );

		// general stats
		int getWidth( void ) { return canvasWidth; }
		int getHeight( void ) { return canvasHeight; }
		int getScale( void ) { return displayScale; }
		void setScale( int newScale );
		void getName( wxString& name );
		void adjustCanvasWidth( void );
		void adjustCanvasHeight( void );

		// some drawing-related convenience methods
		void draw( wxDC& dc );
		void redraw( void );
		void blit( wxDC& dc );
		void blit( void );

		// pens and fonts for drawing different items
		//holds info about a pen
		PenInfo PenArray[numPens];
	private:
		void PutInPenArray(int index,gmtkPen * penptr, wxString nameabr, wxString namelong)
		{
			PenArray[index].pen = penptr;
			PenArray[index].nameabr = nameabr;
			PenArray[index].namelong = namelong;
		}
	public:
		gmtkPen detPen;
		gmtkPen randPen;
		gmtkPen switchingPen;
		gmtkPen det_randPen;
		gmtkPen det_switchPen;
		gmtkPen rand_switchPen;
		gmtkPen det_rand_switchPen;
		gmtkPen highlightPen;
		gmtkPen frameBorderPen;
		gmtkPen chunkBorderPen;
		gmtkPen controlPointPen;
		gmtkPen nodePen;
		gmtkPen gridPen;
		gmtkPen boundingBoxPen;

		wxFont labelFont;

		// Are we drawing ... ?
		bool getViewCPs( void ) { return drawCPs; }
		bool getViewLines( void ) { return drawLines; }
		bool getViewSplines( void ) { return drawSplines; }
		bool getViewArrowHeads( void ) { return drawArrowHeads; }
		bool getViewNodes( void ) { return drawNodes; }
		bool getViewGrids( void ) { return drawGrids; }
		bool getViewDirectLines( void ) { return drawDirectLines; }
		bool getViewFrameSeps( void ) { return drawFrameSeps; }
		bool getViewNodeNames( void ) { return drawNodeNames; }
		bool getViewFrameNames( void ) { return drawFrameNames; }
		bool getViewToolTips( void ) { return drawToolTips; }
		bool getViewBoundingBox( void ) { return drawBoundingBox; }
		bool getViewSelectBox( void ) { return drawSelectBox; }

		// toggle drawing ...
		void toggleViewCPs( void );
		void toggleViewLines( void );
		void toggleViewSplines( void );
		void toggleViewArrowHeads( void );
		void toggleViewNodes( void );
		void toggleViewGrids( void );
		void toggleViewDirectLines( void );
		void toggleViewFrameSeps( void );
		void toggleViewNodeNames( void );
		void toggleViewFrameNames( void );
		void toggleViewToolTips( void );
		void toggleViewBoundingBox( void );
		void toggleViewSelectBox( void );

		void hideSelectedLabels( void );
		void showAllNodeLabels( void );
		void showAllFrameLabels( void );
		void showAllNodes( void );

		// Ask the user what font they want to use.
		void changeFont( void );

		bool getSnapToGrid( void ) { return snapToGrid; }
		void toggleSnapToGrid( void );

		//DECLARE_DYNAMIC_CLASS(StructPage)
		DECLARE_EVENT_TABLE()
	private:
			struct StructPage_state {
				PenInfo PenArray[numPens];
				wxFont labelFont;
				// What do we draw?
				bool drawSelectBox;
				bool drawCPs;
				bool drawLines;
				bool drawSplines;
				bool drawArrowHeads;
				bool drawNodes;
				bool drawGrids;
				bool drawDirectLines;
				bool drawFrameSeps;
				bool drawNodeNames;
				bool drawFrameNames;
				bool drawToolTips;
				bool drawBoundingBox;
				// How big?
				int displayScale;
				long canvasWidth;
				long canvasHeight;

				bool snapToGrid;

				std::vector< VizNode* > nodes; // node positions
				std::vector< NameTag* > nodeNameTags; // node nametags
				std::vector< std::vector< VizArc* > > arcs; // arcs[a][b] is the
				// information about the arc/spline from node a to node b (or NULL)
				std::vector< VizSep* > frameEnds; // positions of dividers between frames
				std::vector< NameTag* > frameNameTags; // frame nametags
			};
	public:
			void save_undo();
			bool undo();
			bool redo();
	private:
			void save_state(std::stack<StructPage_state *> &);
			bool restore_state(std::stack<StructPage_state *> &);
			//this contains states we've been in
			std::stack<StructPage_state *> undo_stack;
			std::stack<StructPage_state *> redo_stack;
			wxFrame *parentFrame;
			wxNotebook *parentNotebook;
			wxBitmap *content; // the drawing buffer
			std::map< RVInfo::rvParent, unsigned int > nameVizNodeMap;
			wxRect selectBox;

			// What do we draw?
			bool drawSelectBox;
			bool drawCPs;
			bool drawLines;
			bool drawSplines;
			bool drawArrowHeads;
			bool drawNodes;
			bool drawGrids;
			bool drawDirectLines;
			bool drawFrameSeps;
			bool drawNodeNames;
			bool drawFrameNames;
			bool drawToolTips;
			bool drawBoundingBox;

			// How big?
			int displayScale;
			long canvasWidth;
			long canvasHeight;
			long rightMostItemX( void );
			long bottomMostItemY( void );

			bool snapToGrid;

			void initNodes( void );
			void initArcs( void );

			// utility methods
			VizArc* newArc(int i, int j);
			VizArc* findArcOwning( ControlPoint *cp, int& index );
			bool deleteSelectedCps( void );
			bool crossesNode( wxCoord x0, wxCoord x1 );
			bool itemMovingInBounds(int frame, wxCoord x_cur, wxCoord x_new);
			bool inBounds( wxCoord x, wxCoord y );
			int nodeUnderPt(wxPoint pt);

			// things related to the structure (str) file
			wxString strFile;
			FileParser *fp;

			// things related to the position (gvp) file
			wxString gvpFile;
			wxString gvpFile_old;	//the previous gvpfile, used when saving as
			/// needs saving?
			bool gvpDirty;
			/// gvp file had peculiarities?
			bool gvpAborted;
			// reads the gvp info into the config map
			void fillMap( void );
			map<const wxString, wxString, ltstr> config;
			std::vector< VizNode* > nodes; // node positions
			std::vector< NameTag* > nodeNameTags; // node nametags
			std::vector< std::vector< VizArc* > > arcs; // arcs[a][b] is the
			// information about the arc/spline from node a to node b (or NULL)
			std::vector< VizSep* > frameEnds; // positions of dividers between frames
			std::vector< NameTag* > frameNameTags; // frame nametags
			std::vector< int > firstNodeInFrame; // nodes in each frame;
			int firstChunkFrame;
			int firstEpilogueFrame;
			int numFrames;
			bool data_parsed; //indicates if the str file parsed or not
};


/// Represents anything that can be selected
class Selectable {
	public:
		// Is it currently selected?
		virtual bool getSelected( void ) { return selected; }
		// Tell it whether or not to be selected.
		virtual void setSelected( bool newSelected ) { selected = visible ? newSelected : false;}
		virtual void toggleSelected( void ) { setSelected(!getSelected()); }
		// If it can be selected, it must be somewhere. These tell you
		// if a click landed on it...
		virtual bool onMe( const wxPoint& pt ) = 0;
		// or if it's in a given rectangle (e.g. the selection rectangle)
		virtual bool inRect( const wxRect& rect ) = 0;
		// constructor: items start out unselected
		Selectable() { selected = false; visible = true;}
		// destructor: no dnamically allocated memory to free
		virtual 	~Selectable() {  }

		virtual void setVisible(bool newVisible){visible = newVisible;}
		void toggleVisible(){setVisible(!visible);}
		bool isVisible(){return visible;}

	protected:
		bool selected;
		bool visible;
};


/// Represents a node's nametag
class NameTag : public Selectable {
	private:
		int frame;
	public:
		/// absolute position (not relative to node)
		wxPoint pos;
		/// automatically updated size
		wxPoint size;
		/// the text displayed
		wxString name;
		// constructor
		NameTag( const wxPoint& newPos, const wxString& newName );
		//copy constructor
		NameTag( const NameTag &);
		// drawing
		void draw( wxDC *dc );
		// selectable method overrides
		virtual bool onMe( const wxPoint& pt );
		virtual bool inRect( const wxRect& rect );
		int GetFrame() { return frame;}
		void SetFrame(int frameNum) { frame = frameNum;}
};

/// Represents a node 
class VizNode : public Selectable {
	public:
		/// absolute position of the center of the circle
		wxPoint center;
		/// a copy of the RVInfo entry from the FileParser
		RVInfo *rvi;
		/// the parent
		StructPage *page;
		/// name and frame number
		RVInfo::rvParent rvId;
		/// the node's label
		NameTag nametag;
		/// pointer to a generic window from a failed attempt at tooltips
		wxWindow *tipWin;
		// constructor
		VizNode( const wxPoint& newPos, RVInfo *newRvi, StructPage *parentPage );
		//copy constructor
		VizNode( const VizNode &);
		// destructor (noop for now)
		~VizNode( void );
		// draw the node
		void draw( wxDC *dc );
		// get frame number
		int GetFrame(){return rvId.second;}
		// Selectable methods
		virtual void setSelected( bool newSelected );
		virtual bool onMe( const wxPoint& pt );
		virtual bool inRect( const wxRect& rect ) { return rect.Inside(center); }

		//highlight functions
		bool isHighlighted() {return (highlight_state != off);}
		void setHighlightState(highlight_states state) {highlight_state = state;}
		highlight_states getHighlightState() {return highlight_state;}
		highlight_states getNextHighlightState();

		//overload set Visible
		void setVisible( bool newVisible );
		
	private:
		highlight_states highlight_state;
};


/// Represents a control point
class ControlPoint : public Selectable {
	public:
		/// absolute position
		wxPoint pos;
		/// the parent arc owning this control point
		VizArc *arc;
		// constructor
		ControlPoint( const wxPoint& pt );
		// selectable methods
		virtual bool onMe( const wxPoint& pt );
		virtual bool inRect( const wxRect& rect ) { return rect.Inside(pos); }
};

/// Represents an arc
class VizArc {
	public:
		/// for some wxDC methods
		wxList *points;
		/// much more user-friendly and with all the necessary info
		std::vector< ControlPoint* > *cps;
		/// the parent
		StructPage *page;

		//these methods set the type of this Arc
		void setRand(void){comb_type |= random;}
		void setSwitching(void){comb_type |= switching;}
		void setDeterministic(void){comb_type |= determin;}

		// constructor
		VizArc( std::vector< ControlPoint* > *newCps, StructPage *newPage );
		// copy constructor
		VizArc(const VizArc&);
		// destructor (to delete the wxList of points)
		~VizArc( void );
		enum { DRAW_ARCS = 1<<0, DRAW_CPS = 1<<1 };
		// draw itself
		void draw( wxDC *dc, int drawFlags );
		//visiblity
		void setVisible(bool newVisible){visible = newVisible;}
		bool isVisible(){return visible;}
		//highlight state
		bool isHighlighted() {return highlighted;}
		void setHighlighted(bool newHighlighted) {highlighted = newHighlighted;}
		void toggleHighlighted() {setHighlighted(!highlighted);}
	private:
		// is it visible
		bool visible;
		bool highlighted;
		// combined type {deterministic, random, switching, deterministic/random,
		// deterministic/switching....}
		enum type {undefined, determin = 1, random, det_rand, switching, det_switch, rand_switch,det_rand_switch};
		unsigned int comb_type;
};


/// Represents a separator between frames
class VizSep : public Selectable {
public:
	/// only need an absolute x position for this
	wxCoord x;
	/// Does it separate the chunk from the prologue or epilogue?
	enum FrameSepType { PROLOGUE, BEGIN_CHUNK, CHUNK, END_CHUNK, EPILOGUE };
	FrameSepType sepType;
	/// the parent
	StructPage *page;
	// constructor
	VizSep( wxCoord xNew, StructPage *newPage, FrameSepType newSepType );
	//copy constructor
	VizSep(const VizSep &);
	// draw itself
	void draw( wxDC *dc );
	// selectable methods
	virtual bool onMe( const wxPoint& pt );
	virtual bool inRect( const wxRect& rect ) { return false; }
};

/// This is what wxWidgets uses for printing, previewing, etc.
class GmtkPrintout : public wxPrintout {
public:
	// constructor
	GmtkPrintout( StructPage *newPage,
		  const wxChar *title = wxT("gmtkViz printout") );
	// What to do when a page needs to be printed.
	bool OnPrintPage(int page);
	// We only have one page
	bool HasPage(int pageNum) { return (pageNum == 1); }
	// Tell the system about our pages.
	void GetPageInfo( int *minPage, int *maxPage,
			  int *selPageFrom, int *selPageTo );
	// actually set the scales and draw
	void DrawPageOne(wxDC *dc);
private:
	/// The StructPage whose draw() method we use.
	StructPage *page;
};

//this class provides a simple help window for gmtkViz
class GmtkHelp : public wxFrame {
	public:
		GmtkHelp(wxWindow* parent, wxWindowID id) : wxFrame(parent, id, "GmktViz Help", 
				wxDefaultPosition, wxDefaultSize, wxDEFAULT_FRAME_STYLE,
				"gmtkVizHelp") {}
		void DisplayContents(void){Show(true); Raise();}
		void doLayout();
		DECLARE_EVENT_TABLE()
	private:
		//we need 2 events because the close button and window close button
		//create different types of events
		void OnClose(wxCloseEvent &event){ Show(false); }
		void OnCloseButton(wxCommandEvent &event){ Show(false); }
};

//event table for GmtkHelp
BEGIN_EVENT_TABLE(GmtkHelp, wxFrame)
	EVT_CLOSE(GmtkHelp::OnClose)
	EVT_BUTTON(wxID_CLOSE, GmtkHelp::OnCloseButton)
END_EVENT_TABLE()

/// This is the main window
class GFrame: public wxFrame {
public:
	/// all the menu event IDs
	enum {
		MENU_FILE_NEW = 1006,
		MENU_FILE_OPEN,
		MENU_FILE_SAVE,
		MENU_FILE_SAVEAS,
		MENU_FILE_PAGESETUP,
		MENU_FILE_PRINT,
		MENU_FILE_PRINT_EPS,
		MENU_FILE_CLOSE,
		MENU_FILE_EXIT,
	MENU_EDIT_SELECTALL,
	MENU_EDIT_SNAPTOGRID,
	MENU_EDIT_CANVASWIDTH,
	MENU_EDIT_CANVASHEIGHT,
	MENU_EDIT_COPYFRAMELAYOUT,
	MENU_EDIT_COPYPARTITIONLAYOUT,
	MENU_EDIT_UNDO,
	MENU_EDIT_REDO,
	MENU_HELP,
	MENU_VIEW_HIDELABELS,
	MENU_VIEW_SHOW_NODE_LABELS,
	MENU_VIEW_SHOW_FRAME_LABELS,
	MENU_VIEW_ALL_NODES_VISIBLE,
	MENU_VIEW_CPS,
	MENU_VIEW_LINES,
	MENU_VIEW_SPLINES,
	MENU_VIEW_ARROW_HEADS,
	MENU_VIEW_NODES,
	MENU_VIEW_GRIDS,
	MENU_VIEW_DIRECT_LINES,
	MENU_VIEW_FRAME_SEPS,
	MENU_VIEW_NODE_NAMES,
	MENU_VIEW_FRAME_NAMES,
	MENU_VIEW_TOOLTIPS,
	MENU_VIEW_BOUNDING_BOX,
	MENU_ZOOM_BEGIN,
	//took these out because they caused crashes with dotted lines
//	MENU_ZOOM_2_pow_neg_4dot000,

//	MENU_ZOOM_2_pow_neg_3dot000,
//	MENU_ZOOM_2_pow_neg_2dot000,
	MENU_ZOOM_2_pow_neg_1dot000,
	MENU_ZOOM_2_pow_neg_0dot500,
	MENU_ZOOM_2_pow_neg_0dot250,
	MENU_ZOOM_2_pow_neg_0dot125,
	MENU_ZOOM_2_pow_pos_0dot000,
	MENU_ZOOM_2_pow_pos_0dot125,
	MENU_ZOOM_2_pow_pos_0dot250,
	MENU_ZOOM_2_pow_pos_0dot500,
	MENU_ZOOM_2_pow_pos_1dot000,
	MENU_ZOOM_2_pow_pos_2dot000,
	//took these out because they cause a loss of responce with large
	//graphics
//	MENU_ZOOM_2_pow_pos_3dot000,
//	MENU_ZOOM_2_pow_pos_4dot000,
	MENU_ZOOM_END,
	MENU_CUSTOMIZE_FONT,
	MENU_CUSTOMIZE_PENS
	/*
	MENU_CUSTOMIZE_PENS_BEGIN,
	MENU_CUSTOMIZE_SWITCHING_PEN,
	MENU_CUSTOMIZE_SWITCHING_PEN_COLOR,
	MENU_CUSTOMIZE_SWITCHING_PEN_WIDTH,
	MENU_CUSTOMIZE_SWITCHING_PEN_STYLE,
	MENU_CUSTOMIZE_CONDITIONAL_PEN,
	MENU_CUSTOMIZE_CONDITIONAL_PEN_COLOR,
	MENU_CUSTOMIZE_CONDITIONAL_PEN_WIDTH,
	MENU_CUSTOMIZE_CONDITIONAL_PEN_STYLE,
	MENU_CUSTOMIZE_BOTH_PEN,
	MENU_CUSTOMIZE_BOTH_PEN_COLOR,
	MENU_CUSTOMIZE_BOTH_PEN_WIDTH,
	MENU_CUSTOMIZE_BOTH_PEN_STYLE,
	MENU_CUSTOMIZE_HIGHLIGHT_PEN,
	MENU_CUSTOMIZE_HIGHLIGHT_PEN_COLOR,
	MENU_CUSTOMIZE_HIGHLIGHT_PEN_WIDTH,
	MENU_CUSTOMIZE_HIGHLIGHT_PEN_STYLE,
	MENU_CUSTOMIZE_FRAMEBORDER_PEN,
	MENU_CUSTOMIZE_FRAMEBORDER_PEN_COLOR,
	MENU_CUSTOMIZE_FRAMEBORDER_PEN_WIDTH,
	MENU_CUSTOMIZE_FRAMEBORDER_PEN_STYLE,
	MENU_CUSTOMIZE_CHUNKBORDER_PEN,
	MENU_CUSTOMIZE_CHUNKBORDER_PEN_COLOR,
	MENU_CUSTOMIZE_CHUNKBORDER_PEN_WIDTH,
	MENU_CUSTOMIZE_CHUNKBORDER_PEN_STYLE,
	MENU_CUSTOMIZE_CONTROLPOINT_PEN,
	MENU_CUSTOMIZE_CONTROLPOINT_PEN_COLOR,
	MENU_CUSTOMIZE_CONTROLPOINT_PEN_WIDTH,
	MENU_CUSTOMIZE_CONTROLPOINT_PEN_STYLE,
	MENU_CUSTOMIZE_NODE_PEN,
	MENU_CUSTOMIZE_NODE_PEN_COLOR,
	MENU_CUSTOMIZE_NODE_PEN_WIDTH,
	MENU_CUSTOMIZE_NODE_PEN_STYLE,
	MENU_CUSTOMIZE_GRID_PEN,
	MENU_CUSTOMIZE_GRID_PEN_COLOR,
	MENU_CUSTOMIZE_GRID_PEN_WIDTH,
	MENU_CUSTOMIZE_GRID_PEN_STYLE,
	MENU_CUSTOMIZE_BOUNDING_BOX_PEN,
	MENU_CUSTOMIZE_BOUNDING_BOX_PEN_COLOR,
	MENU_CUSTOMIZE_BOUNDING_BOX_PEN_WIDTH,
	MENU_CUSTOMIZE_BOUNDING_BOX_PEN_STYLE,
	MENU_CUSTOMIZE_PENS_END
	*/
	};

	/**
	 * Constructor. Creates a new GFrame.
	 */
	GFrame( wxWindow* parent, int id, const wxString& title,
		const wxPoint& pos=wxDefaultPosition,
		const wxSize& size=wxDefaultSize,
		long style=wxDEFAULT_FRAME_STYLE );

	void file(wxString &fileName, bool gvpFormat);

	/**
	 * Processes menu File->New
	 */
	void OnMenuFileNew(wxCommandEvent &event);

	/**
	 * Processes menu File->Open
	 */
	void OnMenuFileOpen(wxCommandEvent &event);

	/**
	 * Processes menu File->Save
	 */
	void OnMenuFileSave(wxCommandEvent &event);

	/**
	 * Processes menu File->Save As
	 */
	void OnMenuFileSaveas(wxCommandEvent &event);

	/**
	 * Processes menu File->Print
	 */
	void OnMenuFilePrint(wxCommandEvent &event);
	
	/**
	 * Processes menu File->PrintEPS
	 */
	void OnMenuFilePrintEPS(wxCommandEvent &event);

	/**
	 * Processes menu File->Page Setup
	 */
	void OnMenuFilePageSetup(wxCommandEvent &event);

	/**
	 * Processes menu File->Close
	 */
	void OnMenuFileClose(wxCommandEvent &event);

	/**
	 * Processes menu File->Exit
	 */
	void OnMenuFileExit(wxCommandEvent &event);

	/**
	 * Processes close window events
	 */
	void OnClose(wxCloseEvent &event);

	void OnMenuEditSelectAll(wxCommandEvent &event);
	void OnMenuEditSnaptogrid(wxCommandEvent &event);
	void OnMenuEditCanvaswidth(wxCommandEvent &event);
	void OnMenuEditCanvasheight(wxCommandEvent &event);
	void OnMenuEditCopyframelayout(wxCommandEvent &event);
	void OnMenuEditCopypartitionlayout(wxCommandEvent &event);
	void OnMenuEditUndo(wxCommandEvent &event);
	void OnMenuEditRedo(wxCommandEvent &event);

	void OnMenuHelp(wxCommandEvent &event);

	// Handle events from the View menu to toggle drawing various items
	void OnMenuViewHideLabels(wxCommandEvent &event);
	void OnMenuViewShowNodeLabels(wxCommandEvent &event);
	void OnMenuViewShowFrameLabels(wxCommandEvent &event);
	void OnMenuViewAllNodesVisible(wxCommandEvent &event);
	void OnMenuViewCPs(wxCommandEvent &event);
	void OnMenuViewLines(wxCommandEvent &event);
	void OnMenuViewSplines(wxCommandEvent &event);
	void OnMenuViewArrowHeads(wxCommandEvent &event);
	void OnMenuViewNodes(wxCommandEvent &event);
	void OnMenuViewGrids(wxCommandEvent &event);
	void OnMenuViewDirectLines(wxCommandEvent &event);
	void OnMenuViewFrameSeps(wxCommandEvent &event);
	void OnMenuViewNodeNames(wxCommandEvent &event);
	void OnMenuViewFrameNames(wxCommandEvent &event);
	void OnMenuViewToolTips(wxCommandEvent &event);
	void OnMenuViewBoundingBox(wxCommandEvent &event);

	// Handle events from the Zoom menu to change the scale/zoom
	void OnMenuZoom(wxCommandEvent &event);

	// Handle events from the Customize menu to alter how items are drawn
	void OnMenuCustomizeFont(wxCommandEvent &event);
	void OnMenuCustomizePen(wxCommandEvent &event);

	// Do this when a different notebook page is chosen
	void OnNotebookPageChanged(wxNotebookEvent &event);

	//update the menu to reflect the given structpage
	void UpdateMenuChecks(StructPage *);

private:
	// initialize the status bar
	void set_properties();
	// arrange and insert the about page
	void do_layout();
	// places to keep the print and page setup settings
	wxPrintData printData;
	wxPageSetupData pageSetupData;

	GmtkHelp * helpWindow;

protected:
	// the widgets associated with this GFrame
	wxMenuBar* MainVizWindow_menubar;
	wxStatusBar* MainVizWindow_statusbar;
	wxStaticText* about_label;
	wxTextCtrl* about_info;
	wxPanel* about_pane;
	wxNotebook* struct_notebook;

	DECLARE_EVENT_TABLE()
};

//global
GFrame* MainVizWindow;

/// represents the app as a whole to the wxWidgets system
class GMTKStructVizApp: public wxApp {
public:
	// kind of like main(), except this just gets things started
	bool OnInit();
};

/*** implementation starts here ***/

IMPLEMENT_APP(GMTKStructVizApp)

#include "arguments.h"
bool help = false;
bool print_version_and_exit = false;
#define MAX_OBJECTS (5)
char *gvpFileNames[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL};
char *strFileNames[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL};
char *cppCommandOptions = NULL;
int verbosity;

Arg Arg::Args[] = {
	Arg( "cppCommandOptions", Arg::Opt, cppCommandOptions,
	 "Additional CPP command line" ),
	Arg( "gvpFile", Arg::Opt, gvpFileNames, "position file",
		Arg::ARRAY, MAX_OBJECTS ),
	Arg( "strFile", Arg::Opt, strFileNames, "structure file",
		Arg::ARRAY, MAX_OBJECTS ),
	Arg("verbosity",Arg::Opt,verbosity,"Verbosity (0 <= v <= 100) of informational/debugging msgs"),
	Arg( "help", Arg::Tog, help, "print this message" ),
	Arg()
};
/**
 *******************************************************************
 * Called implicitly by the wxWidgets system, this method shows
 * the main window and gets the the event loop started.
 *
 * \pre This method should not be called explicitly.
 *
 * \post The main window is shown and the program is running.
 *
 * \note Shows the main window and generally starts the program.
 *
 * \return true
 *******************************************************************/
bool GMTKStructVizApp::OnInit()
{
	bool parse_was_ok = Arg::parse(argc, argv);
	if(help || !parse_was_ok) {
		Arg::usage();
		return false;
	}

	//init the global Style Maps
	//these are used to map the gvp file's text representation of style
	//to wxWidgets numerical representation and visa versa
	styleStringToIntMap["wxSOLID"] = wxSOLID;
   styleIntToStringMap[wxSOLID] = "wxSOLID";
	styleToStyleChoicesIndexMap[wxSOLID] = 0;
	styleStringToIntMap["wxTRANSPARENT"] = wxTRANSPARENT;
   styleIntToStringMap[wxTRANSPARENT] = "wxTRANSPARENT";
	styleToStyleChoicesIndexMap[wxTRANSPARENT] = 1;
	styleStringToIntMap["wxDOT"] = wxDOT;
   styleIntToStringMap[wxDOT] = "wxDOT";
	styleToStyleChoicesIndexMap[wxDOT] = 2;
	styleStringToIntMap["wxLONG_DASH"] = wxLONG_DASH;
   styleIntToStringMap[wxLONG_DASH] = "wxLONG_DASH";
	styleToStyleChoicesIndexMap[wxLONG_DASH] = 3;
	styleStringToIntMap["wxSHORT_DASH"] = wxSHORT_DASH;
   styleIntToStringMap[wxSHORT_DASH] = "wxSHORT_DASH";
	styleToStyleChoicesIndexMap[wxSHORT_DASH] = 4;
	styleStringToIntMap["wxDOT_DASH"] = wxDOT_DASH;
   styleIntToStringMap[wxDOT_DASH] = "wxDOT_DASH";
	styleToStyleChoicesIndexMap[wxDOT_DASH] = 5;
	styleStringToIntMap["wxSTIPPLE"] = wxSTIPPLE;
   styleIntToStringMap[wxSTIPPLE] = "wxSTIPPLE";
	styleToStyleChoicesIndexMap[wxSTIPPLE] = 6;
	styleStringToIntMap["wxUSER_DASH"] = wxUSER_DASH;
   styleIntToStringMap[wxUSER_DASH] = "wxUSER_DASH";
	styleToStyleChoicesIndexMap[wxUSER_DASH] = 7;
	styleStringToIntMap["wxBDIAGONAL_HATCH"] = wxBDIAGONAL_HATCH;
   styleIntToStringMap[wxBDIAGONAL_HATCH] = "wxBDIAGONAL_HATCH";
	styleToStyleChoicesIndexMap[wxBDIAGONAL_HATCH] = 8;
	styleStringToIntMap["wxCROSSDIAG_HATCH"] = wxCROSSDIAG_HATCH;
   styleIntToStringMap[wxCROSSDIAG_HATCH] = "wxCROSSDIAG_HATCH";
	styleToStyleChoicesIndexMap[wxCROSSDIAG_HATCH] = 9;
	styleStringToIntMap["wxFDIAGONAL_HATCH"] = wxFDIAGONAL_HATCH;
   styleIntToStringMap[wxFDIAGONAL_HATCH] = "wxFDIAGONAL_HATCH";
	styleToStyleChoicesIndexMap[wxFDIAGONAL_HATCH] = 10;
	styleStringToIntMap["wxCROSS_HATCH"] = wxCROSS_HATCH;
   styleIntToStringMap[wxCROSS_HATCH] = "wxCROSS_HATCH";
	styleToStyleChoicesIndexMap[wxCROSS_HATCH] = 11;
	styleStringToIntMap["wxHORIZONTAL_HATCH"] = wxHORIZONTAL_HATCH;
   styleIntToStringMap[wxHORIZONTAL_HATCH] = "wxHORIZONTAL_HATCH";
	styleToStyleChoicesIndexMap[wxHORIZONTAL_HATCH] = 12;
	styleStringToIntMap["wxVERTICAL_HATCH"] = wxVERTICAL_HATCH;
   styleIntToStringMap[wxVERTICAL_HATCH] = "wxVERTICAL_HATCH";
	styleToStyleChoicesIndexMap[wxVERTICAL_HATCH] = 13;

	//Set up the Pen Defaults (default within the program)
	
	//detPen
	PenDefault_t * temp_default = new PenDefault_t;
	temp_default->color = *wxBLACK;
	temp_default->style = wxSOLID;
	temp_default->width = ACTUAL_SCALE;
	PenDefaultMap["detPen"] = temp_default;
	//randPen
	temp_default = new PenDefault_t;
	temp_default->color = *wxRED;
	temp_default->style = wxSOLID;
	temp_default->width = ACTUAL_SCALE;
	PenDefaultMap["randPen"] = temp_default;
	//switchingPen
	temp_default = new PenDefault_t;
	temp_default->color = *wxBLUE;
	temp_default->style = wxDOT;
	temp_default->width = ACTUAL_SCALE;
	PenDefaultMap["switchingPen"] = temp_default;
	//det_randPen
	temp_default = new PenDefault_t;
	temp_default->color = *wxGREEN;
	temp_default->style = wxSOLID;
	temp_default->width = ACTUAL_SCALE;
	PenDefaultMap["det_randPen"] = temp_default;
	//det_switchPen
	temp_default = new PenDefault_t;
	temp_default->color = *wxGREEN;
	temp_default->style = wxSOLID;
	temp_default->width = ACTUAL_SCALE;
	PenDefaultMap["det_switchPen"] = temp_default;
	//rand_switchPen
	temp_default = new PenDefault_t;
	temp_default->color = *wxGREEN;
	temp_default->style = wxSOLID;
	temp_default->width = ACTUAL_SCALE;
	PenDefaultMap["rand_switchPen"] = temp_default;
	//det_rand_switchPen
	temp_default = new PenDefault_t;
	temp_default->color = *wxGREEN;
	temp_default->style = wxSOLID;
	temp_default->width = ACTUAL_SCALE;
	PenDefaultMap["det_rand_switchPen"] = temp_default;
	//highlightPen
	temp_default = new PenDefault_t;
	temp_default->color = wxColour(255,0,255);
	temp_default->style = wxSOLID;
	temp_default->width = ACTUAL_SCALE + 2;
	PenDefaultMap["highlightPen"] = temp_default;
	//frameBorderPen
	temp_default = new PenDefault_t;
	temp_default->color = *wxLIGHT_GREY;
	temp_default->style = wxDOT;
#ifdef __WXMSW__
	temp_default->width = 1;
#else
	temp_default->width = ACTUAL_SCALE;
#endif
	PenDefaultMap["frameBorderPen"] = temp_default;
	//chunkBorderPen
	temp_default = new PenDefault_t;
	temp_default->color = *wxBLACK;
	temp_default->style = wxSOLID;
	temp_default->width = ACTUAL_SCALE;
	PenDefaultMap["chunkBorderPen"] = temp_default;
	//controlPointPen
	temp_default = new PenDefault_t;
	temp_default->color = *wxRED;
	temp_default->style = wxSOLID;
	temp_default->width = ACTUAL_SCALE + 1;
	PenDefaultMap["controlPointPen"] = temp_default;
	//nodePen
	temp_default = new PenDefault_t;
	temp_default->color = *wxBLACK;
	temp_default->style = wxSOLID;
	temp_default->width = ACTUAL_SCALE;
	PenDefaultMap["nodePen"] = temp_default;
	//gridPen
	temp_default = new PenDefault_t;
	temp_default->color = *wxLIGHT_GREY;
	temp_default->style = wxSOLID;
	temp_default->width = ACTUAL_SCALE;
	PenDefaultMap["gridPen"] = temp_default;
	//boundingBoxPen
	temp_default = new PenDefault_t;
	temp_default->color = *wxBLACK;
	temp_default->style = wxSOLID;
	temp_default->width = ACTUAL_SCALE;
	PenDefaultMap["boundingBoxPen"] = temp_default;

	//set the default font
	defaultFont = wxFont(12*ACTUAL_SCALE, wxMODERN, wxNORMAL,wxNORMAL);

	//load up any data that is in ~/.Xdefaults
	//these are user defined defaults, they override the program
	//defaults, they are not saved to the gvp file
#ifndef __WXMSW__
	Display *display_name;
	display_name = XOpenDisplay(NULL);

	//the interator for the PenDefaultMap
	map<const wxString, PenDefault_t *, ltstr>::iterator pen_iterator;
	wxString PenWidth;
	wxString PenStyle;
	wxString PenColor;
	wxString value;


	//for each default pen see if there is data in ~/.Xdefaults
	for (pen_iterator = PenDefaultMap.begin(); pen_iterator != PenDefaultMap.end(); pen_iterator++){
		//width
		PenWidth = pen_iterator->first;
		PenWidth = PenWidth.Append("Width");
		value = XGetDefault(display_name,"gmtkViz",PenWidth.c_str());
		if (value != wxEmptyString) {
			long w;
			if (value.ToLong(&w)) {
				pen_iterator->second->width = w;
			} else {
				cout << "\"gmtkViz." << PenWidth << " : " << value << "\" is not a valid number (this is probably in ~/.Xdefaults)\n";
			}
		}

		//color
		PenColor = pen_iterator->first;
		PenColor = PenColor.Append("Color");
		value = XGetDefault(display_name,"gmtkViz",PenColor.c_str());
		if (value != wxEmptyString) {
			int r,g,b;
			if(sscanf(value.c_str(),"%i,%i,%i",&r,&g,&b) == 3){
				pen_iterator->second->color = wxColour(r,g,b);
			} else {
				wxColour temp_color = wxColour(value);
				if (temp_color.Ok())
					pen_iterator->second->color = temp_color;
				else {
					cout << "\"gmtkViz." << PenColor << " : " << value << 
						"\" is not a valid color r,g,b or string (ie \"red\") (this is probably in ~/.Xdefaults)\n";
				}
			}
		}

		//style
		PenStyle = pen_iterator->first;
		PenStyle = PenStyle.Append("Style");
		value = XGetDefault(display_name,"gmtkViz",PenStyle.c_str());
		if (value != wxEmptyString) {
			//if there is anything with that Style in the map
			if (styleStringToIntMap.count(value)) {
				pen_iterator->second->style = styleStringToIntMap[value];
			} else {
				cout << "\"gmtkViz." << PenStyle << " : " << value << 
					"\" is not a valid style (like wxSOLID, wxDOT, etc) (this is probably in ~/.Xdefaults)\n";
			}
		}
	}

	//look for a user defined default font

	//font name
	//this will only work for 2.6.1 and above, 2.4.2 doesn't let you set
	//SetNativeFontInfo based on a string
#if wxCHECK_VERSION(2,6,1)
	value = XGetDefault(display_name,"gmtkViz","fontName");
	if (value != wxEmptyString) {
		defaultFont.SetNativeFontInfo(value);
	}
#endif
	//font size
	value = XGetDefault(display_name,"gmtkViz","fontSize");
	if (value != wxEmptyString) {
		long w;
		if (value.ToLong(&w) && w > 0) {
			defaultFont.SetPointSize(w);
		} else {
			cout << "\"gmtkViz.fontSize: " << value << 
				"\" is not a valid number (this is probably in ~/.Xdefaults)\n";
		}
	}
	
	//close the display
   XCloseDisplay(display_name);
#endif
	
	(void) IM::setGlbMsgLevel(verbosity);
	GM_Parms.setMsgLevel(verbosity);

	wxInitAllImageHandlers();
	// MainVizWindow has no parent...
	MainVizWindow = new GFrame( 0, -1,
					wxT("GMTK Structure File Vizualizer"),
					wxDefaultPosition, wxDefaultSize,
					wxDEFAULT_FRAME_STYLE );
	// ...because it's the top level window.
	SetTopWindow(MainVizWindow);
	// Once we show the window, the program is event driven.
	MainVizWindow->Show();
	for (int i = 0; i < MAX_OBJECTS && gvpFileNames[i]; i++) {
		wxString fileName;
		fileName = gvpFileNames[i];
		MainVizWindow->file(fileName, true);
	}
	for (int i = 0; i < MAX_OBJECTS && strFileNames[i]; i++) {
		wxString fileName;
		fileName = strFileNames[i];
		MainVizWindow->file(fileName, false);
	}
	return true;
}

/**
 *******************************************************************
 * Constructs the main window, setting up menus, tabs, event
 * handlers, etc.
 *
 * \param parent A pointer to the parent of this GFrame. The parent
 *	  will delete this when it is deleted.
 * \param id An integer you may use to identify this window. If you
 *	  don't care, you may specify -1 and wxWidgets will assign it an
 *	  internal identifier.
 * \param argParser A parser set up for strFile, gvpFile, and "" whose
 *	  Parse() method has already been successfully called, so that
		we may simply use the Found() and GetParamCount() methods.
 * \param title What would you like to appear in the window's title
 *	  bar?
 * \param pos Where do you want the window to appear?
 * \param size How big should the window be?
 * \param style An OR'd list of flags regarding the appearance of this
 *	  window. See the wxWidgets documentation for more details.
 *
 * \pre Requires valid parent pointer, id (-1 if you don't care), and
 *	  title. The rest can be omitted. Can't be called explicitly
 *	  (because it's a constructor).
 *
 * \post A GFrame exists and is ready to be shown, etc.
 *
 * \note The GFrame gets registered with the parent to be destroyed
 *	  when the parent is closed.
 *
 * \return Nothing.
 *******************************************************************/
GFrame::GFrame( wxWindow* parent, int id, const wxString& title,
		const wxPoint& pos, const wxSize& size, long style )
	: wxFrame( parent, id, title, pos, size, style )
{
	// create widgets
	struct_notebook = new wxNotebook(this, -1, wxDefaultPosition, wxDefaultSize, wxCLIP_CHILDREN);
	about_pane = new wxPanel(struct_notebook, -1);
	MainVizWindow_menubar = new wxMenuBar();

	// A bunch of menu bar stuff
	SetMenuBar(MainVizWindow_menubar);
	// The File menu
	wxMenu* menu_file = new wxMenu();
	menu_file->Append(MENU_FILE_NEW, wxT("&New...\tCtrl+N"), wxT("Create a new placement file (requires an existing structure file)"), wxITEM_NORMAL);
	menu_file->Append(MENU_FILE_OPEN, wxT("&Open...\tCtrl+O"), wxT("Open an existing placement file"), wxITEM_NORMAL);
	menu_file->AppendSeparator();
	menu_file->Append(MENU_FILE_SAVE, wxT("&Save\tCtrl+S"), wxT("Save the current placement file"), wxITEM_NORMAL);
	menu_file->Append(MENU_FILE_SAVEAS, wxT("Save &As...\tCtrl+Shift+S"), wxT("Save the current placement file with a different name"), wxITEM_NORMAL);
	menu_file->AppendSeparator();
	menu_file->Append(MENU_FILE_PAGESETUP, wxT("Page Setup..."), wxT("Set up page size/orientation"), wxITEM_NORMAL);
	menu_file->Append(MENU_FILE_PRINT, wxT("&Print...\tCtrl+P"), wxT("Preview and print the current graph"), wxITEM_NORMAL);
	menu_file->Append(MENU_FILE_PRINT_EPS, wxT("&Print to EPS file...\tCtrl+E"), wxT("Print an Encapsulated PostScript file of the current graph"), wxITEM_NORMAL);
	menu_file->AppendSeparator();
	menu_file->Append(MENU_FILE_CLOSE, wxT("&Close\tCtrl+W"), wxT("Close current placement file"), wxITEM_NORMAL);
	menu_file->Append(MENU_FILE_EXIT, wxT("E&xit\tCtrl+Q"), wxT("Close all files and exit"), wxITEM_NORMAL);
	MainVizWindow_menubar->Append(menu_file, wxT("&File"));
	// These don't make sense until a document is open.
	MainVizWindow_menubar->Enable(MENU_FILE_SAVE, false);
	MainVizWindow_menubar->Enable(MENU_FILE_SAVEAS, false);
	MainVizWindow_menubar->Enable(MENU_FILE_PRINT, false);
	MainVizWindow_menubar->Enable(MENU_FILE_PRINT_EPS, false);
	MainVizWindow_menubar->Enable(MENU_FILE_CLOSE, false);
	// The Edit menu
	wxMenu* menu_edit = new wxMenu();
	menu_edit->Append(MENU_EDIT_UNDO, wxT("Undo"), wxT("Undo the Last Action"), wxITEM_NORMAL);
	menu_edit->Append(MENU_EDIT_REDO, wxT("Redo"), wxT("Redo the Last Action"), wxITEM_NORMAL);
	menu_edit->Append(MENU_EDIT_SELECTALL, wxT("Select All\tCtrl+A"), wxT("Select Everything that is Selectable"), wxITEM_NORMAL);
	menu_edit->Append(MENU_EDIT_SNAPTOGRID, wxT("Snap To Grid"), wxT("Toggle whether items snap to the grids when they are moved"), wxITEM_CHECK);
	menu_edit->Append(MENU_EDIT_CANVASWIDTH, wxT("Canvas Width..."), wxT("Adjust the width of the canvas"), wxITEM_NORMAL);
	menu_edit->Append(MENU_EDIT_CANVASHEIGHT, wxT("Canvas Height..."), wxT("Adjust the height of the canvas"), wxITEM_NORMAL);
	menu_edit->Append(MENU_EDIT_COPYFRAMELAYOUT, wxT("Copy Frame Layout...\tCtrl+F"), wxT("Copy the layout from one frame to another"), wxITEM_NORMAL);
	menu_edit->Append(MENU_EDIT_COPYPARTITIONLAYOUT, wxT("Copy Partition Layout...\tCtrl+G"), wxT("Copy the layout from one partition to another"), wxITEM_NORMAL);
	MainVizWindow_menubar->Append(menu_edit, wxT("Edit"));
	// The View menu
	wxMenu* menu_view = new wxMenu();
	menu_view->Append(MENU_VIEW_HIDELABELS, wxT("Hide Selected Labels"), wxT("Turn off drawing for all currently selected labels"), wxITEM_NORMAL);
	menu_view->Append(MENU_VIEW_SHOW_NODE_LABELS, wxT("Show All Node Labels"), wxT("Turn on drawing for all node labels"), wxITEM_NORMAL);
	menu_view->Append(MENU_VIEW_SHOW_FRAME_LABELS, wxT("Show All Frame Labels"), wxT("Turn on drawing for all frame labels"), wxITEM_NORMAL);
	menu_view->Append(MENU_VIEW_ALL_NODES_VISIBLE, wxT("Make All Nodes Visible\tCtrl+D"), wxT("Make all nodes Visible"), wxITEM_NORMAL);
	menu_view->AppendSeparator();
	menu_view->Append(MENU_VIEW_CPS, wxT("Draw Control Points\t`"), wxT("Toggle display of arc spline control points"), wxITEM_CHECK);
	menu_view->Append(MENU_VIEW_LINES, wxT("Draw Arc Lines\t1"), wxT("Toggle display of straight lines between control points in arcs"), wxITEM_CHECK);
	menu_view->Append(MENU_VIEW_SPLINES, wxT("Draw Arc Splines\t2"), wxT("Toggle display of arc splines"), wxITEM_CHECK);
	menu_view->Append(MENU_VIEW_ARROW_HEADS, wxT("Draw Arrow Heads\t3"), wxT("Toggle display of arrow heads on arcs"), wxITEM_CHECK);
	menu_view->Append(MENU_VIEW_NODES, wxT("Draw Nodes\t4"), wxT("Toggle display of nodes"), wxITEM_CHECK);
	menu_view->Append(MENU_VIEW_GRIDS, wxT("Draw Grids\t5"), wxT("Toggle display of grids"), wxITEM_CHECK);
	menu_view->Append(MENU_VIEW_DIRECT_LINES, wxT("Draw Direct Lines\t6"), wxT("Toggle display of direct straight lines for arcs"), wxITEM_CHECK);
	menu_view->Append(MENU_VIEW_FRAME_SEPS, wxT("Draw Frame Separators\t7"), wxT("Toggle display of frame separators"), wxITEM_CHECK);
	menu_view->Append(MENU_VIEW_NODE_NAMES, wxT("Draw Node Names\t8"), wxT("Toggle display of node names"), wxITEM_CHECK);
	menu_view->Append(MENU_VIEW_FRAME_NAMES, wxT("Draw Frame Names\t9"), wxT("Toggle display of frame names"), wxITEM_CHECK);
	menu_view->Append(MENU_VIEW_BOUNDING_BOX, wxT("Draw Bounding Box\t0"), wxT("Toggle display of bounding box"), wxITEM_CHECK);
	// XXX: menu_view->Append(MENU_VIEW_TOOLTIPS, wxT("Draw Tool Tips"), wxT("Toggle display of tool tips for node names"), wxITEM_CHECK);
	MainVizWindow_menubar->Append(menu_view, wxT("View"));
	// Doesn't make sense unless a document is active, so disable it for now.
	MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("View")), false);
	// The Zoom menu
	wxMenu* menu_zoom = new wxMenu();
	for (int i = 0; i < MENU_ZOOM_END - MENU_ZOOM_BEGIN - 1; i++) {
		wxString zoomStr;
		zoomStr.sprintf("%7.3f", gZoomMap[i]*(ACTUAL_SCALE));
		menu_zoom->Append( i + MENU_ZOOM_BEGIN + 1, zoomStr,
					wxEmptyString, wxITEM_RADIO );
	}
	MainVizWindow_menubar->Append(menu_zoom, wxT("Zoom"));
	// Also doesn't make sense without a document
	MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("Zoom")), false);
	// The Customize menu
	wxMenu* menu_customize = new wxMenu();
	menu_customize->Append( MENU_CUSTOMIZE_FONT, wxT("Change Font..."),
				wxEmptyString, wxITEM_NORMAL );
	menu_customize->Append( MENU_CUSTOMIZE_PENS, wxT("Customize Pens..."),
				"Pop Up the Customize Pens Dialog", wxITEM_NORMAL );
/*
	wxMenu* menu_customize_switching = new wxMenu();
	menu_customize_switching->Append( MENU_CUSTOMIZE_SWITCHING_PEN_COLOR,
					  wxT("Change Color..."),
					  wxT( "Change the color of arcs "
					   "from switching parents" ),
					  wxITEM_NORMAL );
	menu_customize_switching->Append( MENU_CUSTOMIZE_SWITCHING_PEN_WIDTH,
					  wxT("Change Width..."),
					  wxT( "Change the width of arcs "
					   "from switching parents" ),
					  wxITEM_NORMAL );
	menu_customize_switching->Append( MENU_CUSTOMIZE_SWITCHING_PEN_STYLE,
					  wxT("Change Style..."),
					  wxT( "Change the style of arcs "
					   "from switching parents" ),
					  wxITEM_NORMAL );
	menu_customize->Append( MENU_CUSTOMIZE_SWITCHING_PEN,
				wxT("Switching Pen"), menu_customize_switching );
	wxMenu* menu_customize_conditional = new wxMenu();
	menu_customize_conditional->Append( MENU_CUSTOMIZE_CONDITIONAL_PEN_COLOR,
					  wxT("Change Color..."),
					  wxT( "Change the color of arcs "
					   "from conditional parents" ),
					  wxITEM_NORMAL );
	menu_customize_conditional->Append( MENU_CUSTOMIZE_CONDITIONAL_PEN_WIDTH,
					  wxT("Change Width..."),
					  wxT( "Change the width of arcs "
					   "from conditional parents" ),
					  wxITEM_NORMAL );
	menu_customize_conditional->Append( MENU_CUSTOMIZE_CONDITIONAL_PEN_STYLE,
					  wxT("Change Style..."),
					  wxT( "Change the style of arcs "
					   "from conditional parents" ),
					  wxITEM_NORMAL );
	menu_customize->Append( MENU_CUSTOMIZE_CONDITIONAL_PEN,
				wxT("Conditional Pen"),
				menu_customize_conditional );
	wxMenu* menu_customize_both = new wxMenu();
	menu_customize_both->Append( MENU_CUSTOMIZE_BOTH_PEN_COLOR,
					  wxT("Change Color..."),
					  wxT( "Change the color of arcs "
					   "from parents that are both "
					   "switching and conditional" ),
					  wxITEM_NORMAL );
	menu_customize_both->Append( MENU_CUSTOMIZE_BOTH_PEN_WIDTH,
					  wxT("Change Width..."),
					  wxT( "Change the width of arcs "
					   "from parents that are both "
					   "switching and conditional" ),
					  wxITEM_NORMAL );
	menu_customize_both->Append( MENU_CUSTOMIZE_BOTH_PEN_STYLE,
					  wxT("Change Style..."),
					  wxT( "Change the style of arcs "
					   "from parents that are both "
					   "switching and conditional" ),
					  wxITEM_NORMAL );
	menu_customize->Append( MENU_CUSTOMIZE_BOTH_PEN,
				wxT("Both Pen"),
				menu_customize_both );
	wxMenu* menu_customize_highlight = new wxMenu();
	menu_customize_highlight->Append( MENU_CUSTOMIZE_HIGHLIGHT_PEN_COLOR,
					  wxT("Change Color..."),
					  wxT( "Change the color of "
						  "highlighed arcs and nodes"),
					  wxITEM_NORMAL );
	menu_customize_highlight->Append( MENU_CUSTOMIZE_HIGHLIGHT_PEN_WIDTH,
					  wxT("Change Width..."),
					  wxT( "Change the width of "
						  "highlighed arcs and nodes"),
					  wxITEM_NORMAL );
	menu_customize_highlight->Append( MENU_CUSTOMIZE_HIGHLIGHT_PEN_STYLE,
					  wxT("Change Style..."),
					  wxT( "Change the style of "
						  "highlighed arcs and nodes"),
					  wxITEM_NORMAL );
	menu_customize->Append( MENU_CUSTOMIZE_HIGHLIGHT_PEN,
				wxT("Highlight Pen"),
				menu_customize_highlight );
	wxMenu* menu_customize_frameborder = new wxMenu();
	menu_customize_frameborder->Append( MENU_CUSTOMIZE_FRAMEBORDER_PEN_COLOR,
					  wxT("Change Color..."),
					  wxT("Change the color of frameborders"),
					  wxITEM_NORMAL );
	menu_customize_frameborder->Append( MENU_CUSTOMIZE_FRAMEBORDER_PEN_WIDTH,
					  wxT("Change Width..."),
					  wxT("Change the width of frameborders"),
					  wxITEM_NORMAL );
	menu_customize_frameborder->Append( MENU_CUSTOMIZE_FRAMEBORDER_PEN_STYLE,
					  wxT("Change Style..."),
					  wxT( "Change the style of frameborders"),
					  wxITEM_NORMAL );
	menu_customize->Append( MENU_CUSTOMIZE_FRAMEBORDER_PEN,
				wxT("Frameborder Pen"),
				menu_customize_frameborder );
	wxMenu* menu_customize_chunkborder = new wxMenu();
	menu_customize_chunkborder->Append( MENU_CUSTOMIZE_CHUNKBORDER_PEN_COLOR,
					  wxT("Change Color..."),
					  wxT("Change the color of chunkborders" ),
					  wxITEM_NORMAL );
	menu_customize_chunkborder->Append( MENU_CUSTOMIZE_CHUNKBORDER_PEN_WIDTH,
					  wxT("Change Width..."),
					  wxT("Change the width of chunkborders"),
					  wxITEM_NORMAL );
	menu_customize_chunkborder->Append( MENU_CUSTOMIZE_CHUNKBORDER_PEN_STYLE,
					  wxT("Change Style..."),
					  wxT("Change the style of chunkborders"),
					  wxITEM_NORMAL );
	menu_customize->Append( MENU_CUSTOMIZE_CHUNKBORDER_PEN,
				wxT("Chunkborder Pen"),
				menu_customize_chunkborder );
	wxMenu* menu_customize_controlpoint = new wxMenu();
	menu_customize_controlpoint->Append( MENU_CUSTOMIZE_CONTROLPOINT_PEN_COLOR,
					  wxT("Change Color..."),
					  wxT("Change the color of controlpoints"),
					  wxITEM_NORMAL );
	menu_customize_controlpoint->Append( MENU_CUSTOMIZE_CONTROLPOINT_PEN_WIDTH,
					  wxT("Change Width..."),
					  wxT("Change the width of controlpoints"),
					  wxITEM_NORMAL );
	menu_customize_controlpoint->Append( MENU_CUSTOMIZE_CONTROLPOINT_PEN_STYLE,
					  wxT("Change Style..."),
					  wxT("Change the style of controlpoints"),
					  wxITEM_NORMAL );
	menu_customize->Append( MENU_CUSTOMIZE_CONTROLPOINT_PEN,
				wxT("Controlpoint Pen"),
				menu_customize_controlpoint );
	wxMenu* menu_customize_node = new wxMenu();
	menu_customize_node->Append( MENU_CUSTOMIZE_NODE_PEN_COLOR,
					  wxT("Change Color..."),
					  wxT("Change the color of nodes" ),
					  wxITEM_NORMAL );
	menu_customize_node->Append( MENU_CUSTOMIZE_NODE_PEN_WIDTH,
					  wxT("Change Width..."),
					  wxT("Change the width of nodes"),
					  wxITEM_NORMAL );
	menu_customize_node->Append( MENU_CUSTOMIZE_NODE_PEN_STYLE,
					  wxT("Change Style..."),
					  wxT("Change the style of nodes"),
					  wxITEM_NORMAL );
	menu_customize->Append( MENU_CUSTOMIZE_NODE_PEN,
				wxT("Node Pen"),
				menu_customize_node );
	wxMenu* menu_customize_grid = new wxMenu();
	menu_customize_grid->Append( MENU_CUSTOMIZE_GRID_PEN_COLOR,
					  wxT("Change Color..."),
					  wxT("Change the color of grids" ),
					  wxITEM_NORMAL );
	menu_customize_grid->Append( MENU_CUSTOMIZE_GRID_PEN_WIDTH,
					  wxT("Change Width..."),
					  wxT("Change the width of grids"),
					  wxITEM_NORMAL );
	menu_customize_grid->Append( MENU_CUSTOMIZE_GRID_PEN_STYLE,
					  wxT("Change Style..."),
					  wxT("Change the style of grids"),
					  wxITEM_NORMAL );
	menu_customize->Append( MENU_CUSTOMIZE_GRID_PEN,
				wxT("Grid Pen"),
				menu_customize_grid );
	wxMenu* menu_customize_bounding_box = new wxMenu();
	menu_customize_bounding_box->Append( MENU_CUSTOMIZE_BOUNDING_BOX_PEN_COLOR,
					  wxT("Change Color..."),
					  wxT("Change the color of bounding box" ),
					  wxITEM_NORMAL );
	menu_customize_bounding_box->Append( MENU_CUSTOMIZE_BOUNDING_BOX_PEN_WIDTH,
					  wxT("Change Width..."),
					  wxT("Change the width of bounding box"),
					  wxITEM_NORMAL );
	menu_customize_bounding_box->Append( MENU_CUSTOMIZE_BOUNDING_BOX_PEN_STYLE,
					  wxT("Change Style..."),
					  wxT("Change the style of bounding box"),
					  wxITEM_NORMAL );
	menu_customize->Append( MENU_CUSTOMIZE_BOUNDING_BOX_PEN,
				wxT("Bounding Box Pen"),
				menu_customize_bounding_box );
	*/
	MainVizWindow_menubar->Append(menu_customize, wxT("Customize"));

	// The Help menu
	wxMenu* help_menu = new wxMenu();
	help_menu->Append(MENU_HELP, wxT("Help"), wxT("Pop Up the Help Info Window"), wxITEM_NORMAL);
	MainVizWindow_menubar->Append(help_menu, wxT("Help"));
	
	// Again, needs a document to make sense.
	MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("Customize")), false);

	// The status bar (with 2 fields)
	MainVizWindow_statusbar = CreateStatusBar(2);
	// The text displayed at the top of the About tab
	about_label = new wxStaticText(about_pane, -1, wxT("GMTKStructViz version 0.0.1\n\nThe graph visualizer and organizer for GMTK structure files.\n\nwritten by Evan Dower <evantd@ssli.ee.washington.edu>"), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE);
	// The textarea (and text) at the bottom of the About tab
	about_info = new wxTextCtrl(about_pane, -1, wxT("GMTKStructViz reads GMTK structure files and attempts to display their contents semi-intelligently. Since it is not human, it can only have limited success in this domain. Thus the user is permitted to move nodes and edges to organize the graph in a more logical and visually appealing way than GMTKStructViz's original guess."), wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE|wxTE_READONLY);

	// pretty much justinitialize the status bar
	set_properties();
	// arrange the About tab and add it to the notebook
	do_layout();

	//create the Help Window
	helpWindow = new GmtkHelp(this,-1);
	helpWindow->doLayout();
}

/**
 *******************************************************************
 * Initialize the status bar. I'm thinking of getting rid of this
 * method.
 *
 * \pre Can be called any time after the frame and its status bar have
 *	  been created.
 *
 * \post The status bar text may be altered.
 *
 * \note The status bar text may be altered.
 *
 * \return void
 *******************************************************************/
void GFrame::set_properties()
{
	//SetTitle(wxT("GMTK Structure File Visualizer"));
	/* The second (right-most) field is 16 pixels wide, the first
	 * (left-most) field takes the rest of the space. (Menu tips,
	 * etc. show up in the left-most field. I use the right-most field
	 * to indicate whether the file needs saving or not.) */
	int MainVizWindow_statusbar_widths[] = { -1, 16 };
	MainVizWindow_statusbar->SetStatusWidths(2,MainVizWindow_statusbar_widths);
	const wxString MainVizWindow_statusbar_fields[] = {
	wxT("About GMTKStructViz"),
		wxEmptyString
	};
	for(int i = 0; i < MainVizWindow_statusbar->GetFieldsCount(); ++i) {
		MainVizWindow_statusbar->SetStatusText(MainVizWindow_statusbar_fields[i], i);
	}
}

/**
 *******************************************************************
 * Create and arrange the About tab and add it to the notebook.
 *
 * \pre The frame, notebook, and about pane must exist.
 *
 * \post The about pane is modified and added to the notebook.
 *
 * \note The about pane is modified and added to the notebook.
 *
 * \return void
 *******************************************************************/
void GFrame::do_layout()
{
	wxBoxSizer* MainVizWindow_sizer = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer* about_sizer = new wxBoxSizer(wxVERTICAL);
	wxStaticBoxSizer* about_label_static_sizer = new wxStaticBoxSizer(new wxStaticBox(about_pane, -1, wxT("")), wxHORIZONTAL);
	about_label_static_sizer->Add( about_label, 1,
				   wxALL|wxALIGN_CENTER_HORIZONTAL, 2 );
	about_sizer->Add( about_label_static_sizer, 1,
			  wxEXPAND|wxALIGN_CENTER_HORIZONTAL, 0 );
	about_sizer->Add(about_info, 2, wxALL|wxEXPAND, 2);
	about_pane->SetAutoLayout(true);
	about_pane->SetSizer(about_sizer);
	about_sizer->Fit(about_pane);
	about_sizer->SetSizeHints(about_pane);
	struct_notebook->AddPage(about_pane, wxT("About"));
	MainVizWindow_sizer->Add(new wxNotebookSizer(struct_notebook), 1, wxEXPAND, 0);
	SetAutoLayout(true);
	SetSizer(MainVizWindow_sizer);
	MainVizWindow_sizer->Fit(this);
	MainVizWindow_sizer->SetSizeHints(this);
	Layout();
}

/* *** event handling *** */
// This is where all of GFrame's event handlers are hooked up
BEGIN_EVENT_TABLE(GFrame, wxFrame)
	EVT_MENU(MENU_FILE_NEW, GFrame::OnMenuFileNew)
	EVT_MENU(MENU_FILE_OPEN, GFrame::OnMenuFileOpen)
	EVT_MENU(MENU_FILE_SAVE, GFrame::OnMenuFileSave)
	EVT_MENU(MENU_FILE_SAVEAS, GFrame::OnMenuFileSaveas)
	EVT_MENU(MENU_FILE_PAGESETUP, GFrame::OnMenuFilePageSetup)
	EVT_MENU(MENU_FILE_PRINT, GFrame::OnMenuFilePrint)
	EVT_MENU(MENU_FILE_PRINT_EPS, GFrame::OnMenuFilePrintEPS)
	EVT_MENU(MENU_FILE_CLOSE, GFrame::OnMenuFileClose)
	EVT_MENU(MENU_FILE_EXIT, GFrame::OnMenuFileExit)
	EVT_MENU(MENU_EDIT_SELECTALL, GFrame::OnMenuEditSelectAll)
	EVT_MENU(MENU_EDIT_SNAPTOGRID, GFrame::OnMenuEditSnaptogrid)
	EVT_MENU(MENU_EDIT_CANVASWIDTH, GFrame::OnMenuEditCanvaswidth)
	EVT_MENU(MENU_EDIT_CANVASHEIGHT, GFrame::OnMenuEditCanvasheight)
	EVT_MENU(MENU_EDIT_COPYFRAMELAYOUT, GFrame::OnMenuEditCopyframelayout)
	EVT_MENU(MENU_EDIT_COPYPARTITIONLAYOUT, GFrame::OnMenuEditCopypartitionlayout)
	EVT_MENU(MENU_EDIT_UNDO, GFrame::OnMenuEditUndo)
	EVT_MENU(MENU_EDIT_REDO, GFrame::OnMenuEditRedo)
	EVT_MENU(MENU_HELP, GFrame::OnMenuHelp)
	EVT_MENU(MENU_VIEW_HIDELABELS, GFrame::OnMenuViewHideLabels)
	EVT_MENU(MENU_VIEW_SHOW_NODE_LABELS, GFrame::OnMenuViewShowNodeLabels)
	EVT_MENU(MENU_VIEW_SHOW_FRAME_LABELS, GFrame::OnMenuViewShowFrameLabels)
	EVT_MENU(MENU_VIEW_ALL_NODES_VISIBLE, GFrame::OnMenuViewAllNodesVisible)
	EVT_MENU(MENU_VIEW_CPS, GFrame::OnMenuViewCPs)
	EVT_MENU(MENU_VIEW_LINES, GFrame::OnMenuViewLines)
	EVT_MENU(MENU_VIEW_SPLINES, GFrame::OnMenuViewSplines)
	EVT_MENU(MENU_VIEW_ARROW_HEADS, GFrame::OnMenuViewArrowHeads)
	EVT_MENU(MENU_VIEW_NODES, GFrame::OnMenuViewNodes)
	EVT_MENU(MENU_VIEW_GRIDS, GFrame::OnMenuViewGrids)
	EVT_MENU(MENU_VIEW_DIRECT_LINES, GFrame::OnMenuViewDirectLines)
	EVT_MENU(MENU_VIEW_FRAME_SEPS, GFrame::OnMenuViewFrameSeps)
	EVT_MENU(MENU_VIEW_NODE_NAMES, GFrame::OnMenuViewNodeNames)
	EVT_MENU(MENU_VIEW_FRAME_NAMES, GFrame::OnMenuViewFrameNames)
	EVT_MENU(MENU_VIEW_TOOLTIPS, GFrame::OnMenuViewToolTips)
	EVT_MENU(MENU_VIEW_BOUNDING_BOX, GFrame::OnMenuViewBoundingBox)
	EVT_MENU_RANGE(MENU_ZOOM_BEGIN+1, MENU_ZOOM_END-1, GFrame::OnMenuZoom)
	EVT_MENU(MENU_CUSTOMIZE_FONT, GFrame::OnMenuCustomizeFont)
	EVT_MENU(MENU_CUSTOMIZE_PENS, GFrame::OnMenuCustomizePen)
//	EVT_MENU_RANGE(MENU_CUSTOMIZE_PENS_BEGIN+1, MENU_CUSTOMIZE_PENS_END-1, GFrame::OnMenuCustomizePen)
	EVT_NOTEBOOK_PAGE_CHANGED(wxID_ANY, GFrame::OnNotebookPageChanged)
	EVT_CLOSE(GFrame::OnClose)
END_EVENT_TABLE()


/**
 *******************************************************************
 * Make a new StructPage for the given structure or position.
 *
 * \param fileName The file to use.
 * \param gvpFormat Is it a position file (as opposed to a structure)?
 *
 * \pre The GFrame needs to be fully initialized.
 *
 * \post A new page is added for this file.
 *
 * \note A new page is added for this file.
 *
 * \return void
 *******************************************************************/
void
GFrame::file(wxString &fileName, bool gvpFormat)
{
	// This will add itself to the notebook and be destroyed when
	// the notebook is destroyed.
	StructPage *page = new StructPage( struct_notebook, -1, this,
			struct_notebook, fileName, gvpFormat );
	//see if we've parsed sucessfully
	if (page->parsedSucessfully()){
		int w, h;
		GetSize(&w, &h);
		SetSize( w<page->getWidth()/ACTUAL_SCALE+25 ?
			 page->getWidth()/ACTUAL_SCALE+25 : w,
			 h<page->getHeight()/ACTUAL_SCALE+100 ?
			 page->getHeight()/ACTUAL_SCALE+100 : h );
		// We won't get an event for this new notebook page coming to
		// the front, so we'll just pretend we did
		wxNotebookEvent dummy;
		OnNotebookPageChanged(dummy);
	} else {
		delete page;
	}
}

/**
 *******************************************************************
 * Prompt the user for a structure file and create a StructPage
 * for it. The StructPage will add iteslf to the notebook.
 *
 * \param event Ignored.
 *
 * \pre Program must be fully initialized.
 *
 * \post A new page may have been created and added to the notebook.
 *
 * \note A new page may have been created and added to the notebook.
 *
 * \return void
 *******************************************************************/
void GFrame::OnMenuFileNew(wxCommandEvent &event)
{
	// The position file will be new, but they still have to open a
	// structure file.
	wxFileDialog dlg(this,
			 "Find the desired structure file", "", "",
			 "All files|*"
			 "|Structure Files|*.str"
			 "|Text Files|*.txt;*.text",
			 wxOPEN | wxCHANGE_DIR, wxDefaultPosition);

	dlg.SetFilterIndex(1); // show only .str's by default

	// put the dialog up (as a modal dialog) and only continue if the
	// user clicked OK
	if ( dlg.ShowModal() == wxID_OK ) {
		wxString fileName;
		fileName = dlg.GetPath();
		file(fileName, false);
	}
}

/**
 *******************************************************************
 * Prompt the user for a position file and create a StructPage
 * for it. The StructPage will add iteslf to the notebook.
 *
 * \param event Ignored.
 *
 * \pre Program must be fully initialized.
 *
 * \post A new page may have been created and added to the notebook.
 *
 * \note A new page may have been created and added to the notebook.
 *
 * \return void
 *******************************************************************/
void GFrame::OnMenuFileOpen(wxCommandEvent &event)
{
	wxFileDialog dlg(this,
			 "Find the desired position file", "", "",
			 "All files|*"
			 "|Position Files|*.gvp"
			 "|Text Files|*.txt;*.text",
			 wxOPEN | wxCHANGE_DIR, wxDefaultPosition);

	dlg.SetFilterIndex(1); // show only .gvp's by default

	// Show the dialog modally and only continue if the user dismissed
	// the dialog by clicking OK
	if ( dlg.ShowModal() == wxID_OK ) {
		wxString fileName;
		fileName = dlg.GetPath();
		file(fileName, true);
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage to do the actual
 * saving and user prompting.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post The file may have been saved.
 *
 * \note The file may have been saved.
 *
 * \return void
 *******************************************************************/
void GFrame::OnMenuFileSave(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage)
		curPage->Save();
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage to do the actual
 * saving and user prompting.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post The file may have been saved.
 *
 * \note The file may have been saved.
 *
 * \return void
 *******************************************************************/
void GFrame::OnMenuFileSaveas(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage)
		curPage->SaveAs();
}

/**
 *******************************************************************
 * Display a page setup dialog and alter the page setup data.
 *
 * \param event Ignored.
 *
 * \pre Program must be fully initialized.
 *
 * \post The print and page setup data may have been altered.
 *
 * \note The print and page setup data may have been altered.
 *
 * \return void
 *******************************************************************/
void GFrame::OnMenuFilePageSetup(wxCommandEvent &event)
{
	pageSetupData = printData;

	wxPageSetupDialog pageSetupDialog(this, &pageSetupData);

	if (pageSetupDialog.ShowModal() == wxID_OK) {
		printData = pageSetupDialog.GetPageSetupData().GetPrintData();
		pageSetupData = pageSetupDialog.GetPageSetupData();
	}
}

/**
 *******************************************************************
 * Pass the buck to a GmtkPrintout which will then pass the buck
 * to the appropriate StructPage to do the actual
 * printing, previewing, and user prompting.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post The file may have been printed and some printing preferences
 *	  may have been altered..
 *
 * \note The file may have been printed and some printing preferences
 *	  may have been altered..
 *
 * \return void
 *******************************************************************/


//This class subclasses wxPrintPrview and disables zoom (for the time being)
class nozoomwxPrintPreview : public wxPrintPreview
{
	public:
		void SetZoom(int zoom);
		nozoomwxPrintPreview(GmtkPrintout*, GmtkPrintout*, wxPrintDialogData* = NULL);
};

nozoomwxPrintPreview::nozoomwxPrintPreview(GmtkPrintout* printout, 
		GmtkPrintout* printoutForPrinting, wxPrintDialogData* data) 
: wxPrintPreview(printout, printoutForPrinting, data)
{
	wxPrintPreview::SetZoom(200);
};

void nozoomwxPrintPreview::SetZoom(int zoom)
{
	//do nothing
	new wxTipWindow(NULL, "Zooming has been known to cause crashes when\n"
			"non-solid lines are used so zooming has been disabled.", 500, NULL);
	return;
};

void GFrame::OnMenuFilePrint(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage) {
		// get a title for the printout
		wxString name;
		curPage->getName(name);
			
		/* Pass two printout objects: for preview, and possible
		 * printing. The printout object is what does much of the
		 * work. It in turn calls the StructPage's draw() method to do
		 * the actual drawing. */
		wxPrintDialogData printDialogData(printData);

//if we're using 2.6.1 or greater we don't get crashes on zoom
//so there is no need to disable it
#if wxCHECK_VERSION(2,6,1)
		wxPrintPreview *preview =
			new wxPrintPreview(new GmtkPrintout(curPage, name),
						new GmtkPrintout(curPage, name),
						&printDialogData);
#else
		nozoomwxPrintPreview *preview =
			new nozoomwxPrintPreview(new GmtkPrintout(curPage, name),
						new GmtkPrintout(curPage, name),
						&printDialogData);
#endif
		if (!preview->Ok()) {
			delete preview;
			wxMessageBox(_T("There was a problem previewing.\n"
						"Perhaps your current printer is not set correctly?"), 
					wxT("Previewing"), wxOK);
			return;
		}
		// Now that the preview is set up, show it in a preview frame.
		
		//XXX this allows for dotted lines (not ideal)
		//preview->SetZoom(200);

		wxPreviewFrame *frame =
			new wxPreviewFrame(preview, this, wxT("gmtkViz Print Preview"),
						wxPoint(100, 100), wxSize(600, 650));
		frame->Centre(wxBOTH);
		frame->Initialize();
		frame->Show(TRUE);
	}
}

/**
 *******************************************************************
 * Get a file name from the User, print a temporary PostScript file
 * execute "cat <tempPostScriptFile> | ps2eps > <filefromUser>"
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post The file may have been printed 
 *
 * \note The file may have been printed 
 *
 * \return void
 *******************************************************************/

void GFrame::OnMenuFilePrintEPS(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage) {

		wxFileDialog eps_file_dialog(this, "Save to EPS", "", "", "*.eps", wxSAVE | wxOVERWRITE_PROMPT);
		if(eps_file_dialog.ShowModal() == wxID_CANCEL)
			return; //canceled

		//see if the bounding box is being viewed, if not, display it
		//but keep track so we can revert to the previous settings later
		bool had_bounding_box = curPage->getViewBoundingBox();
		if(!had_bounding_box)
			curPage->toggleViewBoundingBox();

		// get a title for the printout
		wxString name;
		curPage->getName(name);
      GmtkPrintout printout(curPage,name);
		wxString temp_file_name = wxGetTempFileName( wxT("gmtkviz") );

		//give the temp file name to the printData
		printData.SetFilename(temp_file_name);
		printData.SetPrintMode(wxPRINT_MODE_FILE);
		wxPrintDialogData printDialogData(printData);
		//set up the printer make it print to file
		printDialogData.SetPrintToFile(true);
		printDialogData.SetPrintData(printData);

      wxPrinter printer(&printDialogData);

//XXX Print direct to file doesn't work in 2.6.1 so we need to pop up the
//print dialog
#if wxCHECK_VERSION(2,6,1)
		wxMessageBox(_T("In order to print to EPS you simply have to press OK, Save, Yes when the print dialog"
					" shows up, don't change any settings"),
				wxT("Printing to EPS"), wxOK);
		if(!printer.Print(this, &printout, true)){
			wxMessageBox(_T("There was a problem printing to an EPS file.\n"
						"gmtkViz having trouble printing to file."),
					wxT("Error Printing to EPS"), wxOK);
			return;
		}
#else
		if(!printer.Print(this, &printout, false)){
			wxMessageBox(_T("There was a problem printing to an EPS file.\n"
						"gmtkViz having trouble printing to file."),
					wxT("Error Printing to EPS"), wxOK);
			return;
		}
#endif
		
		//this is the command being sent to the shell (single quote the file name)
		string ps2eps_cmd = "cat ";
		ps2eps_cmd.append(printer.GetPrintDialogData().GetPrintData().GetFilename());
		ps2eps_cmd.append(" | ps2eps > '");

		//quote single quotes to be safe
		string temp_path = eps_file_dialog.GetPath().c_str();
		size_t match_loc = 0;
		while(string::npos != (match_loc = temp_path.find("'", match_loc))){
			temp_path.replace(match_loc,1,"'\"'\"'"); //replace ' with '"'"'
			match_loc += 5;	//skip what we put in
		}
		ps2eps_cmd.append(temp_path);
		ps2eps_cmd.append("'");
		
		//execute the command
		if(system(ps2eps_cmd.c_str()) != 0){
			//print error if needed
			wxMessageBox(_T("There was a problem printing to an EPS file.\n"
						"Perhaps ps2eps is not in your path?\n"
						"Check the shell output for more information."), 
					wxT("Error Printing to EPS"), wxOK);
		}

		//get rid of temp file
		wxRemoveFile(temp_file_name);
		
		//revert to previous settings
		if(!had_bounding_box)
			curPage->toggleViewBoundingBox();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage to do any
 * necessary prompting and saving. Delete the page if it says we
 * can.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post The file may have been saved and/or closed.
 *
 * \note The file may have been saved and/or closed.
 *
 * \return void
 *******************************************************************/
void GFrame::OnMenuFileClose(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage) {
		// delete the page only if it says it's okay
		if (curPage->RequestClose()) {
			struct_notebook->DeletePage(curPageNum);
		}
	}
}

/**
 *******************************************************************
 * Just attempt to close the frame.
 *
 * \param event Ignored.
 *
 * \pre Program should probably be fully initialized.
 *
 * \post Files may have been saved and/or closed and the program may
 *	  have exited.
 *
 * \note Files may have been saved and/or closed and the program may
 *	  have exited.
 *
 * \return void
 *******************************************************************/
void GFrame::OnMenuFileExit(wxCommandEvent &event)
{
	// Just try to close the window. The details will be handled in
	// the event handler OnClose().
	Close(false);
}

/**
 *******************************************************************
 * Loop through the notebook pages, requesting that any
 * StructPages close themselves. Pass the buck to the appropriate
 * StructPages to do any necessary prompting and saving. Delete
 * the pages when they say we can, and cancel if any say we can't.
 *
 * \param event Ignored.
 *
 * \pre The program should be fully initialized. This is meant to be
 *	  called automatically by wxWidgets when an attempt is made to
 *	  close the window.
 *
 * \post Files may have been saved and/or closed and the main window
 *	  may have been closed.
 *
 * \note Files may have been saved and/or closed and the main window
 *	  may have been closed.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnClose(wxCloseEvent& event)
{
	bool destroy = true;

	// There are some cases where we're not allowed to veto the event.
	if ( event.CanVeto() ) {
		StructPage *curPage;
		// for each tab (going backward since we delete them)
		for ( int pageNum = struct_notebook->GetPageCount() - 1;
				pageNum >= 0 && destroy;
				pageNum-- ) {
			curPage = dynamic_cast<StructPage*>
				(struct_notebook->GetPage(pageNum));
			/* If it couldn't be casted to a StructPage, then curPage
			 * will be NULL. */
			if (curPage) {
				// only delete the tab if it says we can
				if (curPage->RequestClose()) {
					struct_notebook->DeletePage(pageNum);
				} else {
					// otherwise the user cancelled so we abort
					event.Veto();
					destroy = false;
				}
			}
		}
	}
	// if we didn't abort (weren't cancelled), then destroy the frame
	if ( destroy ) {
		Destroy();
	}
}
/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * select all of its selectable items.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then all of its selectable items
 * will be selected
 *
 * \note If a StructPage was in front, then all of its selectable items
 * will be selected
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuEditSelectAll(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->setAllSelected(true);
		curPage->redraw();
		curPage->blit();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle snapping to the grids.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  snaps to the grids.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  snaps to the grids.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuEditSnaptogrid(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->toggleSnapToGrid();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * prompt the user for a new width.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then its width may have changed.
 *
 * \note If a StructPage was in front, then its width may have changed.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuEditCanvaswidth(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->adjustCanvasWidth();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * prompt the user for a new height.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then its height may have changed.
 *
 * \note If a StructPage was in front, then its height may have changed.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuEditCanvasheight(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->adjustCanvasHeight();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * copy frame layout
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has asked the user
 *   about which frame to copy to which and porbably moved nodes around.
 *
 * \note If a StructPage was in front, then it has asked the user
 *   about which frame to copy to which and porbably moved nodes around.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuEditCopyframelayout(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->copyFrameLayout();
		/* TODO: make frameCopyDialog happen from here
		 * each StructPage should have it's own dialog and if there is one being displayed
		 * when a new page is loaded then the old dialog should be replaced with the new one
		 * We also need to keep track of when a structpage is closed so that if a dialog is
		 * being displayed it will be either closed or replaced by the dialog for the structpage
		 * that'll be shown next
		 *
		if (frameCopyDialog != NULL) {
		//if the page has changed then we need to display a different copyDialog
		} else {
			//do something like this
			frameCopyDialog = curPage->getFrameCopyDialog();
		}
			
		*/
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * copy partition layout
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has asked the user
 *   about which partition to copy to which and porbably moved nodes around.
 *
 * \note If a StructPage was in front, then it has asked the user
 *   about which partition to copy to which and porbably moved nodes around.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuEditCopypartitionlayout(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage)
		curPage->copyPartitionLayout();
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * undo the last action
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front then the last action will be undone
 * 	(if there is one)
 *
 * \note If a StructPage was in front then the last action will be undone
 * 	(if there is one)
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuEditUndo(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		if (!curPage->undo()){
			//make beep sound
			printf("\a");
			fflush(stdout);
		} else {
			curPage->redraw();
			curPage->blit();
		}
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * redo the last action (undo the last undo)
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front then the last action will be redone
 * 	(if there is one)
 *
 * \note If a StructPage was in front then the last action will be redone
 * 	(if there is one)
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuEditRedo(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		if (!curPage->redo()){
			//make beep sound
			printf("\a");
			fflush(stdout);
		} else {
			curPage->redraw();
			curPage->blit();
		}
	}
}


/**
 *******************************************************************
 * Pop Up the Help Window
 *
 * \param event Ignored.
 *
 * \pre The GFrame Should Be fully Initialized.
 *
 * \post A Window will pop Up with Help info.
 *
 * \note
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuHelp(wxCommandEvent &event)
{
	helpWindow->DisplayContents();
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it not to
 * draw the selected labels.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has turned off drawing
 *  for the selected labels.
 *
 * \note If a StructPage was in front, then it has turned off drawing
 *  for the selected labels.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewHideLabels(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->hideSelectedLabels();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * draw all node labels.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has turned on drawing
 *  for all node labels.
 *
 * \note If a StructPage was in front, then it has turned on drawing
 *  for all node labels.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewShowNodeLabels(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->showAllNodeLabels();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * draw all frame labels.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has turned on drawing
 *  for all frame labels.
 *
 * \note If a StructPage was in front, then it has turned on drawing
 *  for all frame labels.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewShowFrameLabels(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->showAllFrameLabels();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * draw all nodes.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has turned on drawing
 *  for all nodes.
 *
 * \note If a StructPage was in front, then it has turned on drawing
 *  for all nodes.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewAllNodesVisible(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->showAllNodes();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle the drawing of control points
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  draws control points and redrawn itself.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  draws control points and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewCPs(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->toggleViewCPs();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle the drawing of straight line segments from control
 * point to control point
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  draws straight line segments and redrawn itself.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  draws straight line segments and redrawn itself.
 *
 * \return void
 */
void
GFrame::OnMenuViewLines(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->toggleViewLines();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle the drawing of arc splines
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  draws arc splines segments and redrawn itself.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  draws arc splines segments and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewSplines(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->toggleViewSplines();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle the drawing of arrow heads
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  draws arrow heads and redrawn itself.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  draws arrow heads and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewArrowHeads(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->toggleViewArrowHeads();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle the drawing of nodes
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  draws nodes and redrawn itself.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  draws nodes and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewNodes(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->toggleViewNodes();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle the drawing of grids
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  draws grids and redrawn itself.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  draws grids and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewGrids(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->toggleViewGrids();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle the drawing of direct lines from node to node
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  draws direct lines and redrawn itself.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  draws direct lines and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewDirectLines(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->toggleViewDirectLines();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle the drawing of frame separators
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  draws frame separators and redrawn itself.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  draws frame separators and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewFrameSeps(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->toggleViewFrameSeps();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle the drawing of node names
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  draws node names and redrawn itself.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  draws node names and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewNodeNames(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage){
		curPage->save_undo();
		curPage->toggleViewNodeNames();
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle the drawing of frame names
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  draws frame names and redrawn itself.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  draws frame names and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewFrameNames(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage)
		curPage->toggleViewFrameNames();
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle the drawing of tooltips. Note that tool tips don't
 * actually work so this doesn't really do anything.
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  draws tooltips and redrawn itself.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  draws tooltips and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewToolTips(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage)
		curPage->toggleViewToolTips();
}
/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * toggle the drawing of the bounding box
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has toggled whether it
 *	  draws the bounding box and redrawn itself.
 *
 * \note If a StructPage was in front, then it has toggled whether it
 *	  draws the bounding box and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuViewBoundingBox(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage)
		curPage->toggleViewBoundingBox();
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * change its scaling factor
 *
 * \param event The event caused by selecting one of the Zoom menu
 *	  items. At the very least event.GetId() must return the right
 *	  thing.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it has adjusted its scaling
 *	  factor and redrawn itself.
 *
 * \note If a StructPage was in front, then it has adjusted its scaling
 *	  factor and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuZoom(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage) {
		curPage->save_undo();
		int id = event.GetId()/*, scale = curPage->getScale()*/;
		// Since the are radio items, only one will be checked at a time.
		//MainVizWindow_menubar->Check(MENU_ZOOM_BEGIN + scale + 1, false);
		MainVizWindow_menubar->Check(id, true);
		curPage->setScale(id - MENU_ZOOM_BEGIN - 1);
	}
}

/**
 *******************************************************************
 * Pass the buck to the appropriate StructPage telling it to
 * change its font
 *
 * \param event Ignored.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it may have gotten a new
 *	  font from the user and redrawn itself.
 *
 * \note If a StructPage was in front, then it may have gotten a new
 *	  font from the user and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuCustomizeFont(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage)
		curPage->changeFont();
}

/**
 *******************************************************************
 * Modify the StructPage's pen and tell it to repaint.
 *
 * \param event The event caused by selecting one of the Customize menu
 *	  items. At the very least event.GetId() must return the right
 *	  thing.
 *
 * \pre A StructPage should be at the front, but precautions are taken
 *	  in case it isn't, so everything should be fine as long as the
 *	  program is fully initialized.
 *
 * \post If a StructPage was in front, then it may have altered one of
 *	  its pens and redrawn itself.
 *
 * \note If a StructPage was in front, then it may have altered one of
 *	  its pens and redrawn itself.
 *
 * \return void
 *******************************************************************/
void
GFrame::OnMenuCustomizePen(wxCommandEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage) {
		curPage->OnCustomizePen();
#if 0
		int id = event.GetId()/*, scale = curPage->getScale()*/;
		wxPen *thePen = NULL;
		switch (id) {
			case MENU_CUSTOMIZE_SWITCHING_PEN_COLOR:
			case MENU_CUSTOMIZE_SWITCHING_PEN_WIDTH:
			case MENU_CUSTOMIZE_SWITCHING_PEN_STYLE:
				thePen = &curPage->switchingPen;
				break;
			case MENU_CUSTOMIZE_CONDITIONAL_PEN_COLOR:
			case MENU_CUSTOMIZE_CONDITIONAL_PEN_WIDTH:
			case MENU_CUSTOMIZE_CONDITIONAL_PEN_STYLE:
				thePen = &curPage->conditionalPen;
				break;
			case MENU_CUSTOMIZE_BOTH_PEN_COLOR:
			case MENU_CUSTOMIZE_BOTH_PEN_WIDTH:
			case MENU_CUSTOMIZE_BOTH_PEN_STYLE:
				thePen = &curPage->bothPen;
				break;
			case MENU_CUSTOMIZE_HIGHLIGHT_PEN_COLOR:
			case MENU_CUSTOMIZE_HIGHLIGHT_PEN_WIDTH:
			case MENU_CUSTOMIZE_HIGHLIGHT_PEN_STYLE:
				thePen = &curPage->highlightPen;
				break;
			case MENU_CUSTOMIZE_FRAMEBORDER_PEN_COLOR:
			case MENU_CUSTOMIZE_FRAMEBORDER_PEN_WIDTH:
			case MENU_CUSTOMIZE_FRAMEBORDER_PEN_STYLE:
				thePen = &curPage->frameBorderPen;
				break;
			case MENU_CUSTOMIZE_CHUNKBORDER_PEN_COLOR:
			case MENU_CUSTOMIZE_CHUNKBORDER_PEN_WIDTH:
			case MENU_CUSTOMIZE_CHUNKBORDER_PEN_STYLE:
				thePen = &curPage->chunkBorderPen;
				break;
			case MENU_CUSTOMIZE_CONTROLPOINT_PEN_COLOR:
			case MENU_CUSTOMIZE_CONTROLPOINT_PEN_WIDTH:
			case MENU_CUSTOMIZE_CONTROLPOINT_PEN_STYLE:
				thePen = &curPage->controlPointPen;
				break;
			case MENU_CUSTOMIZE_NODE_PEN_COLOR:
			case MENU_CUSTOMIZE_NODE_PEN_WIDTH:
			case MENU_CUSTOMIZE_NODE_PEN_STYLE:
				thePen = &curPage->nodePen;
				break;
			case MENU_CUSTOMIZE_GRID_PEN_COLOR:
			case MENU_CUSTOMIZE_GRID_PEN_WIDTH:
			case MENU_CUSTOMIZE_GRID_PEN_STYLE:
				thePen = &curPage->gridPen;
				break;
			case MENU_CUSTOMIZE_BOUNDING_BOX_PEN_COLOR:
			case MENU_CUSTOMIZE_BOUNDING_BOX_PEN_WIDTH:
			case MENU_CUSTOMIZE_BOUNDING_BOX_PEN_STYLE:
				thePen = &curPage->boundingBoxPen;
				break;
			default:
				return;
		}
		wxColour newColor;
		int newWidth;
		int newStyleNum;
		int penStyle;
		static const wxString choices[] = { wxT("Solid"),
							wxT("Transparent"),
							wxT("Dotted"),
							wxT("Dashed (long dashes)"),
							wxT("Dashed (shart dashes)"),
							wxT("Dotted and Dashed") };
		switch (id) {
			case MENU_CUSTOMIZE_SWITCHING_PEN_COLOR:
			case MENU_CUSTOMIZE_CONDITIONAL_PEN_COLOR:
			case MENU_CUSTOMIZE_BOTH_PEN_COLOR:
			case MENU_CUSTOMIZE_HIGHLIGHT_PEN_COLOR:
			case MENU_CUSTOMIZE_FRAMEBORDER_PEN_COLOR:
			case MENU_CUSTOMIZE_CHUNKBORDER_PEN_COLOR:
			case MENU_CUSTOMIZE_CONTROLPOINT_PEN_COLOR:
			case MENU_CUSTOMIZE_NODE_PEN_COLOR:
			case MENU_CUSTOMIZE_GRID_PEN_COLOR:
			case MENU_CUSTOMIZE_BOUNDING_BOX_PEN_COLOR:
				newColor = wxGetColourFromUser(this, thePen->GetColour());
				if (newColor.Ok()){
					curPage->save_undo();
					thePen->SetColour(newColor);
				}
				//else wxLogMessage("User canceled or invalid color");
				break;
			case MENU_CUSTOMIZE_SWITCHING_PEN_WIDTH:
			case MENU_CUSTOMIZE_CONDITIONAL_PEN_WIDTH:
			case MENU_CUSTOMIZE_BOTH_PEN_WIDTH:
			case MENU_CUSTOMIZE_HIGHLIGHT_PEN_WIDTH:
			case MENU_CUSTOMIZE_FRAMEBORDER_PEN_WIDTH:
			case MENU_CUSTOMIZE_CHUNKBORDER_PEN_WIDTH:
			case MENU_CUSTOMIZE_CONTROLPOINT_PEN_WIDTH:
			case MENU_CUSTOMIZE_NODE_PEN_WIDTH:
			case MENU_CUSTOMIZE_GRID_PEN_WIDTH:
			case MENU_CUSTOMIZE_BOUNDING_BOX_PEN_WIDTH:
				penStyle = thePen->GetStyle();
#ifdef __WXMSW__
				if ( penStyle == wxDOT || penStyle == wxLONG_DASH ||
						penStyle == wxSHORT_DASH || penStyle == wxDOT_DASH ||
						penStyle == wxUSER_DASH ) {
					wxLogMessage("Dots and dashes can only be used with pens of width 1, so you can't change the width until you change the style.");
				} else {
					newWidth = wxGetNumberFromUser( wxT("How wide should the pen be?"),
							wxT("Width: "),
							wxT("Change Pen Width"),
							thePen->GetWidth(), 1, 256 );
					if (newWidth > 0){
						curPage->save_undo();
						thePen->SetWidth(newWidth);
					}
				}
#else
				newWidth = wxGetNumberFromUser( wxT("How wide should the pen be?"),
						wxT("Width: "),
						wxT("Change Pen Width"),
						thePen->GetWidth(), 1, 256 );
				if (newWidth > 0){
					curPage->save_undo();
					thePen->SetWidth(newWidth);
				}
#endif
				break;
			case MENU_CUSTOMIZE_SWITCHING_PEN_STYLE:
			case MENU_CUSTOMIZE_CONDITIONAL_PEN_STYLE:
			case MENU_CUSTOMIZE_BOTH_PEN_STYLE:
			case MENU_CUSTOMIZE_HIGHLIGHT_PEN_STYLE:
			case MENU_CUSTOMIZE_FRAMEBORDER_PEN_STYLE:
			case MENU_CUSTOMIZE_CHUNKBORDER_PEN_STYLE:
			case MENU_CUSTOMIZE_CONTROLPOINT_PEN_STYLE:
			case MENU_CUSTOMIZE_NODE_PEN_STYLE:
			case MENU_CUSTOMIZE_GRID_PEN_STYLE:
			case MENU_CUSTOMIZE_BOUNDING_BOX_PEN_STYLE:
				newStyleNum = wxGetSingleChoiceIndex( wxT("What style would you like?"),
						wxT("Change Pen Style"),
						/* the number of pen choices, if set to 2 we remove the
						 * ability to chose dotted and dashed lines */
						6, choices, this );
				if (newStyleNum >= 0) {
					if ( 2 <= newStyleNum && newStyleNum <= 5 ) {
#if !wxCHECK_VERSION(2,6,1)
						if ( wxNO == wxMessageBox("Lines with dots and/or "
									"dashes are known to cause "
									"problems (crashing) in "
									"some situations (including "
									"printing).\nIf you choose "
									"to use dots and dashes, "
									"please remember to save "
									"frequently.\nAre you sure "
									"you want to change to "
									"this style?", "Confirm",
									wxYES_NO, this) ) {
							break;
						}
#endif
					//save the undo
					curPage->save_undo();

#ifdef __WXMSW__
						if ( thePen->GetWidth() != 1 ) {
							thePen->SetWidth(1);
							wxLogMessage( "Dots and dashes require a pen "
									"width of 1.\nThe width has been "
									"altered accordingly." );
						}
#endif
					}
					switch (newStyleNum) {
						case 0:
							thePen->SetStyle(wxSOLID);
							break;
						case 1:
							thePen->SetStyle(wxTRANSPARENT);
							break;
						case 2:
							thePen->SetStyle(wxDOT);
							break;
						case 3:
							thePen->SetStyle(wxLONG_DASH);
							break;
						case 4:
							thePen->SetStyle(wxSHORT_DASH);
							break;
						case 5:
							thePen->SetStyle(wxDOT_DASH);
							break;
							/*case 6:
							  thePen->SetStyle(wxBDIAGONAL_HATCH);
							  break;
							  case 7:
							  thePen->SetStyle(wxCROSSDIAG_HATCH);
							  break;
							  case 8:
							  thePen->SetStyle(wxFDIAGONAL_HATCH);
							  break;
							  case 9:
							  thePen->SetStyle(wxCROSS_HATCH);
							  break;
							  case 10:
							  thePen->SetStyle(wxHORIZONTAL_HATCH);
							  break;
							  case 11:
							  thePen->SetStyle(wxVERTICAL_HATCH);
							  break;*/
						default:
							return;
					}
				}
				break;
			default:
				return;
		}
		curPage->redraw();
		curPage->blit();
#endif
	}
}

/**
 *******************************************************************
 * Update menus, etc. when the page is changed.
 *
 * \param event Ignored.
 *
 * \pre Program should be fully initialized.
 *
 * \post Menus, etc. should be updated.
 *
 * \note Menus, etc. should be updated.
 *
 * \return void
 *******************************************************************/

void GFrame::UpdateMenuChecks(StructPage *curPage)
{
		MainVizWindow_menubar->Check( MENU_VIEW_CPS,
						  curPage->getViewCPs() );
		MainVizWindow_menubar->Check( MENU_VIEW_LINES,
						  curPage->getViewLines() );
		MainVizWindow_menubar->Check( MENU_VIEW_SPLINES,
						  curPage->getViewSplines() );
		MainVizWindow_menubar->Check( MENU_VIEW_ARROW_HEADS,
						  curPage->getViewArrowHeads() );
		MainVizWindow_menubar->Check( MENU_VIEW_NODES,
						  curPage->getViewNodes() );
		MainVizWindow_menubar->Check( MENU_VIEW_GRIDS,
						  curPage->getViewGrids() );
		MainVizWindow_menubar->Check( MENU_VIEW_DIRECT_LINES,
						  curPage->getViewDirectLines() );
		MainVizWindow_menubar->Check( MENU_VIEW_FRAME_SEPS,
						  curPage->getViewFrameSeps() );
		MainVizWindow_menubar->Check( MENU_VIEW_NODE_NAMES,
						  curPage->getViewNodeNames());
		MainVizWindow_menubar->Check( MENU_VIEW_FRAME_NAMES,
						  curPage->getViewFrameNames());
		MainVizWindow_menubar->Check( MENU_VIEW_BOUNDING_BOX,
						  curPage->getViewBoundingBox() );

		MainVizWindow_menubar->Check( curPage->getScale()+MENU_ZOOM_BEGIN+1,
						  true );
}

void
GFrame::OnNotebookPageChanged(wxNotebookEvent &event)
{
	// figure out which page this is for and pass the buck
	int curPageNum = struct_notebook->GetSelection();
	StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
	// If it couldn't be casted to a StructPage, then curPage will be NULL.
	if (curPage) {
		// pass the buck (this will update the status bar)
		curPage->onComeForward();
		// and do a bunch of stuff on our own too
		// These menu items should be enabled when a document is in front.
		MainVizWindow_menubar->Enable(MENU_FILE_SAVE, true);
		MainVizWindow_menubar->Enable(MENU_FILE_SAVEAS, true);
		MainVizWindow_menubar->Enable(MENU_FILE_PRINT, true);
		MainVizWindow_menubar->Enable(MENU_FILE_PRINT_EPS, true);
		MainVizWindow_menubar->Enable(MENU_FILE_CLOSE, true);
		// and this menu
		MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("Edit")), true);
		MainVizWindow_menubar->Check( MENU_EDIT_SNAPTOGRID,
						  curPage->getSnapToGrid() );
		// and this menu
		MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("View")), true);
		// restore the checked status of each item
		UpdateMenuChecks(curPage);

		// tooltips don't work
		// XXX: MainVizWindow_menubar->Check( MENU_VIEW_TOOLTIPS,
		// curPage->getViewToolTips() );
		// Zoom should be enabled for documents
		MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("Zoom")), true);
		// since they are radio items, keeping only one selected is
		// handled automatically
		/*for (int i = 0, scale = curPage->getScale(); i < MENU_ZOOM_END - MENU_ZOOM_BEGIN - 1; i++) {
			MainVizWindow_menubar->Check( i + MENU_ZOOM_BEGIN + 1,
						  i==scale );
		}*/
//		MainVizWindow_menubar->Check( curPage->getScale()+MENU_ZOOM_BEGIN+1,
//						  true );
		// and the Customize menu should be shown for documents as well
		MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("Customize")), true);
	}
	else {
		// otherwise we set the status bar to some default values
		SetStatusText(wxEmptyString, 1);
		SetStatusText(wxT("About GMTKStructViz"), 0);
		// diable menus and menu items that don't apply to the About tab
		MainVizWindow_menubar->Enable(MENU_FILE_SAVE, false);
		MainVizWindow_menubar->Enable(MENU_FILE_SAVEAS, false);
		MainVizWindow_menubar->Enable(MENU_FILE_PRINT, false);
		MainVizWindow_menubar->Enable(MENU_FILE_PRINT_EPS, false);
		MainVizWindow_menubar->Enable(MENU_FILE_CLOSE, false);
		MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("Edit")), false);
		MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("View")), false);
		MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("Zoom")), false);
		MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("Customize")), false);
	}
}

// I don't claim to understand this
//IMPLEMENT_DYNAMIC_CLASS(StructPage, wxScrolledWindow)

// StructPages handle some events too
BEGIN_EVENT_TABLE(StructPage, wxScrolledWindow)
	EVT_PAINT(StructPage::OnPaint)
	EVT_MOUSE_EVENTS(StructPage::OnMouseEvent)
	EVT_CHAR(StructPage::OnChar)
END_EVENT_TABLE()

/**
 *******************************************************************
 * Constructs a structure graph tab, adds it to the
 * parentNotebook, reads a structure file (and possibly a
 * position file), initializes the positions of nodes and
 * arcs, and updates the status bar.
 *
 * \param parent A pointer to this widget's parent widget so that it
 *	  can be informed of its new child.
 * \param id An integer you can use to identify this widget. If you
 *	  don't care, just give -1 and wxWidgets will make up it's own
 *	  internal identifier.
 * \param parentFrame A point to the wxFrame that this is part of so
 *	  that we can update its status bar.
 * \param parentNotebook The notebook into which we should insert
 *	  ourselves.
 * \param file The filename of either a structure file (.str) or a
 *	  position file (.gvp).
 * \param old Are we opening a pre-existing position file? In other
 *	  words, does \p file reference a position file (as opposed to a
 *	  structure file)?
 *
 * \pre All the parameters should be valid.
 *
 * \post A new tab will be created in \p parentNotebook in which the
 *	  user may manipulate the graph.
 *
 * \note A bunch of memory is allocated. A FileParser is created and
 *	  destroyed. Since FileParser uses globals variables, you
 *	  should avoid having a FileParser around when you call this
 *	  method.
 *
 * \return Nothing.
 *******************************************************************/
StructPage::StructPage(wxWindow *parent, wxWindowID id,
			   wxFrame *parentFrame, wxNotebook *parentNotebook,
			   const wxString &file, bool old)
	: wxScrolledWindow( parent, id, wxDefaultPosition, wxDefaultSize,
			wxSUNKEN_BORDER | wxTAB_TRAVERSAL, _T("") ),
	detPen(*wxBLACK_PEN),
	randPen(*wxBLACK_PEN),
	switchingPen(*wxBLACK_PEN),
	det_randPen(*wxBLACK_PEN),
	det_switchPen(*wxBLACK_PEN),
	rand_switchPen(*wxBLACK_PEN),
	det_rand_switchPen(*wxBLACK_PEN),
	highlightPen(wxPen(wxColour(255,0,255),3,wxSOLID)),//highlight pen is pink
	frameBorderPen(*wxLIGHT_GREY_PEN), chunkBorderPen(*wxBLACK_PEN),
	controlPointPen(*wxRED_PEN), nodePen(*wxBLACK_PEN),
	gridPen(*wxLIGHT_GREY_PEN), 
	boundingBoxPen(*wxBLACK_PEN)
	//labelFont(12*ACTUAL_SCALE, wxMODERN, wxNORMAL,wxNORMAL)
{
	nodes.clear();
	nodeNameTags.clear();
	arcs.clear();
	frameEnds.clear();
	frameNameTags.clear();

	labelFont = defaultFont;
	//since we haven't parsed the data yet we mark it false
	data_parsed = false;
	// This is used later to update the status bar
	this->parentFrame = parentFrame;
	// We add ourselves to this
	this->parentNotebook = parentNotebook;
	// The file starts out clean
	gvpDirty = false;
	// No peculiarities have been found in the gvp file yet
	gvpAborted = false;
	canvasWidth = 1;
	canvasHeight = 1;
	/* No need to draw the selection box since the user can't be
	 * selecting anything quite yet */
	drawSelectBox = false;
	// defaults for optionally drawn items
	drawCPs = true;
	drawLines = false;
	drawSplines = true;
	drawArrowHeads = true;
	drawNodes = true;
	drawGrids = false;
	drawDirectLines = false;
	drawFrameSeps = true;
	drawNodeNames = true;
	drawFrameNames = true;
	// Has no real effect since our tooltips don't work
	drawToolTips = true;
	drawBoundingBox = false;

	snapToGrid = false;

	// The drawing area. If we don't make sure it's NULL, we might
	// accidentally delete it later (cause a segfault).
	content = NULL;

	/* load defaults */
	//pens
	detPen.SetStyle( PenDefaultMap["detPen"]->style);
	detPen.SetWidth( PenDefaultMap["detPen"]->width);
	detPen.SetColour(PenDefaultMap["detPen"]->color);
	randPen.SetStyle( PenDefaultMap["randPen"]->style);
	randPen.SetWidth( PenDefaultMap["randPen"]->width);
	randPen.SetColour(PenDefaultMap["randPen"]->color);
	switchingPen.SetStyle( PenDefaultMap["switchingPen"]->style);
	switchingPen.SetWidth( PenDefaultMap["switchingPen"]->width);
	switchingPen.SetColour(PenDefaultMap["switchingPen"]->color);
	det_randPen.SetStyle( PenDefaultMap["det_randPen"]->style);
	det_randPen.SetWidth( PenDefaultMap["det_randPen"]->width);
	det_randPen.SetColour(PenDefaultMap["det_randPen"]->color);
	det_switchPen.SetStyle( PenDefaultMap["det_switchPen"]->style);
	det_switchPen.SetWidth( PenDefaultMap["det_switchPen"]->width);
	det_switchPen.SetColour(PenDefaultMap["det_switchPen"]->color);
	rand_switchPen.SetStyle( PenDefaultMap["rand_switchPen"]->style);
	rand_switchPen.SetWidth( PenDefaultMap["rand_switchPen"]->width);
	rand_switchPen.SetColour(PenDefaultMap["rand_switchPen"]->color);
	det_rand_switchPen.SetStyle( PenDefaultMap["det_rand_switchPen"]->style);
	det_rand_switchPen.SetWidth( PenDefaultMap["det_rand_switchPen"]->width);
	det_rand_switchPen.SetColour(PenDefaultMap["det_rand_switchPen"]->color);

	highlightPen.SetStyle( PenDefaultMap["highlightPen"]->style);
	highlightPen.SetWidth( PenDefaultMap["highlightPen"]->width);
	highlightPen.SetColour(PenDefaultMap["highlightPen"]->color);
	frameBorderPen.SetStyle( PenDefaultMap["frameBorderPen"]->style);
	frameBorderPen.SetWidth( PenDefaultMap["frameBorderPen"]->width);
	frameBorderPen.SetColour(PenDefaultMap["frameBorderPen"]->color);
	chunkBorderPen.SetStyle( PenDefaultMap["chunkBorderPen"]->style);
	chunkBorderPen.SetWidth( PenDefaultMap["chunkBorderPen"]->width);
	chunkBorderPen.SetColour(PenDefaultMap["chunkBorderPen"]->color);
	controlPointPen.SetStyle( PenDefaultMap["controlPointPen"]->style);
	controlPointPen.SetWidth( PenDefaultMap["controlPointPen"]->width);
	controlPointPen.SetColour(PenDefaultMap["controlPointPen"]->color);
	nodePen.SetStyle( PenDefaultMap["nodePen"]->style);
	nodePen.SetWidth( PenDefaultMap["nodePen"]->width);
	nodePen.SetColour(PenDefaultMap["nodePen"]->color);
	gridPen.SetStyle( PenDefaultMap["gridPen"]->style);
	gridPen.SetWidth( PenDefaultMap["gridPen"]->width);
	gridPen.SetColour(PenDefaultMap["gridPen"]->color);
	boundingBoxPen.SetStyle( PenDefaultMap["boundingBoxPen"]->style);
	boundingBoxPen.SetWidth( PenDefaultMap["boundingBoxPen"]->width);
	boundingBoxPen.SetColour(PenDefaultMap["boundingBoxPen"]->color);
	
//	/* Dots and dashes require width to be one. */
//	switchingPen.SetStyle(wxDOT); /* was wxSOLID */
//	conditionalPen.SetStyle(wxSOLID);
//	bothPen.SetStyle(wxLONG_DASH); /* was wxSOLID */
//	highlightPen.SetStyle(wxSOLID); /* was wxSOLID */
//	frameBorderPen.SetStyle(wxDOT); /* was wxSOLID */
//	chunkBorderPen.SetStyle(wxSOLID);
//	// Scale all these lines up.
//	switchingPen.SetWidth(ACTUAL_SCALE);
//	conditionalPen.SetWidth(ACTUAL_SCALE);
//	bothPen.SetWidth(ACTUAL_SCALE);
//	highlightPen.SetWidth(ACTUAL_SCALE + 2);	//make the highlighted pen thicker
//	frameBorderPen.SetWidth(ACTUAL_SCALE);
//	chunkBorderPen.SetWidth(ACTUAL_SCALE);
//	controlPointPen.SetWidth(ACTUAL_SCALE + 1);
	// Make the corners sharp
	controlPointPen.SetJoin(wxJOIN_MITER);
	nodePen.SetWidth(ACTUAL_SCALE);
	gridPen.SetStyle(wxDOT);
	boundingBoxPen.SetStyle(wxSOLID);
	boundingBoxPen.SetJoin(wxJOIN_MITER);

	//save the state of each pen as the default
	detPen.setDefaults();
	randPen.setDefaults();
	switchingPen.setDefaults();
	det_randPen.setDefaults();
	det_switchPen.setDefaults();
	rand_switchPen.setDefaults();
	det_rand_switchPen.setDefaults();
	highlightPen.setDefaults();
	frameBorderPen.setDefaults();
	chunkBorderPen.setDefaults();
	controlPointPen.setDefaults();
	nodePen.setDefaults();
	gridPen.setDefaults();
	boundingBoxPen.setDefaults();

	//Put data in the pen array
	PutInPenArray( 0,&detPen,"detPen","Deterministic Pen");
	PutInPenArray( 1,&randPen,"randPen","Random Pen");
	PutInPenArray( 2,&switchingPen,"switchingPen","Switching Pen");
	PutInPenArray( 3,&det_randPen,"det_randPen","Deterministic/Random Pen");
	PutInPenArray( 4,&det_switchPen,"det_switchPen","Deterministic/Switching Pen");
	PutInPenArray( 5,&rand_switchPen,"rand_switchPen","Random/Switching Pen");
	PutInPenArray( 6,&det_rand_switchPen,"det_rand_switchPen","Deterministic/Random/Switching Pen");
	PutInPenArray( 7,&highlightPen,"highlightPen","Highlight Pen");
	PutInPenArray( 8,&frameBorderPen,"frameBorderPen","Frameborder Pen");
	PutInPenArray( 9,&chunkBorderPen,"chunkBorderPen","Chunkborder Pen");
	PutInPenArray(10,&controlPointPen,"controlPointPen","Control Point Pen");
	PutInPenArray(11,&nodePen,"nodePen","Node Pen");
	PutInPenArray(12,&gridPen,"gridPen","Grid Pen");
	PutInPenArray(13,&boundingBoxPen,"boundingBoxPen","Bounding Box Pen");

	// scroll 10 pixels at a time
	SetScrollRate( 10, 10 );

	// totally useless for now, but may someday come in handy
	wxToolTip::SetDelay(250);
	wxToolTip::Enable(drawToolTips);

//	// mangle the name and set the tab title and status bar
//	int beginning, ending;
//	if ( (beginning = file.rfind('/') + 1) < 1 )
//		beginning = 0;
//	ending = file.rfind(wxT(".str"));
//	if (ending <= beginning)
//		ending = file.rfind(wxT(".gvp"));
//	if (ending <= beginning)
//		ending = file.size();
//	parentNotebook->InsertPage(parentNotebook->GetPageCount(),
//				   this, file.substr(beginning, ending-beginning),
//				   true);

	/* If this was specified as "old" then the file is a gvp file so
	 * we should populate the config map with info from that.
	 * Otherwise, it's the str file and we don't have a gvp file yet. */
	if (old) {
		gvpFile = file;
		fillMap();
	} else {
		strFile = file;
	}

	if (!gvpAborted && strFile.length()) {
		// this indicates if the strFile location contained in the gvp file is
		// valid below we may mark it invalid
		bool strFile_valid = true;
		//these are just char * versions of the filenames
		char * strFile_cstr = (char *)strFile.mb_str(wxConvUTF8);
		char * gvpFile_cstr = (char *)gvpFile.mb_str(wxConvUTF8);
		/* store the strFile name so that we can use it to in the FileParser
		 * this way we only need to have 1 call to the FileParser, even
		 * if we modify the tempFileName
		 */
		wxString tempFileName(strFile);

		if (old){
			/*
			 * here we try to be intelligent about the location of the strfile
			 * since sometimes it isn't so smart in the gvp file
			 */
			int error_num;
			struct stat strFile_stat;

			error_num = stat(strFile, &strFile_stat);
			//so if there is an error the file probably doesn't exist
			//try to find the file
			if (error_num < 0) {
				// here we try to make a filename for the strfile so that
				// it is relative to the location of the gvpFile
				// since dirname and basename modify their arguments we must
				// create temp copies
				wxString tempGVP = wxString::Format("%s",gvpFile_cstr);
				//wxString tempSTR = wxString::Format("%s",strFile_cstr);
				//= dirname(gvpFile) ++ / ++ strFile
				tempFileName = wxString::Format("%s/%s",
						dirname((char *)tempGVP.mb_str(wxConvUTF8)), 
						strFile_cstr);
						//basename((char *)tempSTR.mb_str(wxConvUTF8)));

				if (strFile_cstr[0] == '/'){
					//the strfile is given as an absolute path and there
					//isn't a valid file there
					strFile_valid = false;
				}
				else if(stat(tempFileName.c_str(), &strFile_stat) != 0){
					//we cannot find it
					perror(tempFileName.c_str());
					strFile_valid = false;
				} 
				/* otherwise we've found the location of the strFile, it was given
				 * relative to the location of the gvpFile this is the ideal
				 * location of the strFile relative to the gvpFile
				 * it is stored in tempFileName
				 */
			}
		} 
		
		if (strFile_valid) {

			/*
			 * Here we fork.
			 * We fork to see if we can parse the data, if we are sucessful then we
			 * parse it again, if not (if the child doesn't exit correctly) then we
			 * print an error and exit.
			 *
			 * XXX:
			 * Later this should be replaced with exceptions, this would be quite
			 * simple, all you'd have to keep would be the sucessful return from
			 * fork statment (but wrap it in a try), and the else statment [when
			 * the child doesn't exit sucessfully] would be the catch statement
			 */
#ifdef __WXMSW__
			//XXX does this actually need to happen?
			//we need to get the filename into a format that cygwin can deal with
			//this is a unix path format
			if (tempFileName.Find("cygwin") >= 0)
				tempFileName.Remove(0,tempFileName.Find("cygwin") + 6);
			tempFileName.Replace("\\","/",true);
#endif
			pid_t pid;
			if ( (pid = fork()) < 0) {
				//cannot fork
				wxString msg;
				msg.sprintf(
						"Could not fork so we cannot make sure that we can parse the data\n"
						"without crashing so we will not parse it at all.");
				wxLogMessage(msg);
			} else if (pid == 0) {
				//the child
				// load up the structure file
				fp = new FileParser(tempFileName, cppCommandOptions);
				// parse the file
				fp->parseGraphicalModel();
				// create the rv variable objects
				fp->createRandomVariableGraph();
				// ensure no loops in graph for all possible unrollings.
				fp->ensureValidTemplate();
				//the fp is of no good to the parent, ditch it
				delete fp;
				fp = NULL;
				_exit(0);	//sucess
			} else {
				int status;
				//the parent
				waitpid(pid, &status, 0);
				if(WEXITSTATUS(status) == EXIT_SUCCESS){
					//exited sucessfully
					// load up the structure file
					
					fp = new FileParser(tempFileName, cppCommandOptions);

					// XXX: prevent this from crashing with parse error
					// parse the file
					fp->parseGraphicalModel();
					// create the rv variable objects
					fp->createRandomVariableGraph();
					// ensure no loops in graph for all possible unrollings.
					fp->ensureValidTemplate();

					// XXX: verify that the str and gvp files go together
					// make up initial positions for the nodes
					initNodes();
					if ( frameEnds.size() ) {
						if ( frameEnds[0]->sepType == VizSep::PROLOGUE ||
								frameEnds[0]->sepType == VizSep::BEGIN_CHUNK )
							frameNameTags.push_back(new NameTag(wxPoint(10*ACTUAL_SCALE,
											10*ACTUAL_SCALE), wxT("frame 0 (P)")));
						else
							frameNameTags.push_back(new NameTag(wxPoint(10*ACTUAL_SCALE,
											10*ACTUAL_SCALE), wxT("frame 0 (C)")));
					}
					for (unsigned int i = 0; i < frameEnds.size(); i++) {
						char partition;
						wxString label;
						if ( frameEnds[i]->sepType == VizSep::PROLOGUE )
							partition = 'P';
						else if ( frameEnds[i]->sepType == VizSep::BEGIN_CHUNK ||
								frameEnds[i]->sepType == VizSep::CHUNK )
							partition = 'C';
						else 
							partition = 'E';
						label.sprintf("frame %d (%c)", i+1, partition);
						frameNameTags.push_back(new NameTag( wxPoint( frameEnds[i]->x +
										10*ACTUAL_SCALE,
										10*ACTUAL_SCALE ),
									label ));
					}
					wxString key, value;
					for (unsigned int i = 0; i < frameNameTags.size(); i++) {
						// move the nametag, if requested
						// x
						key.sprintf( "frames[%d].nametag.pos.x", i );
						value = config[key];
						if (value != wxEmptyString) {
							long xPos;
							if (value.ToLong(&xPos)) {
								frameNameTags[i]->pos.x = xPos;
							} else {
								wxString msg;
								msg.sprintf("frames[%d].nametag.pos.x is not a number", i);
								wxLogMessage(msg);
								gvpAborted = true;
							}
						}
						// y
						key.sprintf( "frames[%d].nametag.pos.y", i );
						value = config[key];
						if (value != wxEmptyString) {
							long yPos;
							if (value.ToLong(&yPos)) {
								frameNameTags[i]->pos.y = yPos;
							} else {
								wxString msg;
								msg.sprintf("frames[%d].nametag.pos.y is not a number", i);
								wxLogMessage(msg);
								gvpAborted = true;
							}
						}
						// visibility
						key.sprintf( "frames[%d].nametag.visible", i );
						value = config[key];
						if (value != wxEmptyString) {
							long visible;
							if (value.ToLong(&visible) && (visible == true || visible == false)) {
								frameNameTags[i]->setVisible(visible);
							} else {
								wxString msg;
								msg.sprintf( "frames[%d].nametag.visible is not a number 0 or 1", i);
								wxLogMessage(msg);
								gvpAborted = true;
							}
						}
					}
					initArcs();

					//for each pen read in the customizations from the gvp if they exist
					for(int i = 0; i < numPens; i++){
						//width
						key.sprintf("%s.width", PenArray[i].nameabr.c_str());
						value = config[key];
						if (value != wxEmptyString) {
							long w;
							if (value.ToLong(&w)) {
								PenArray[i].pen->SetWidth(w);
							} else {
								wxString msg;
								msg.sprintf("%s.width is not a number", PenArray[i].nameabr.c_str());
								wxLogMessage(msg);
								gvpAborted = true;
							}
						}
						//colour
						key.sprintf("%s.colour", PenArray[i].nameabr.c_str());
						value = config[key];
						if (value != wxEmptyString) {
							int r,g,b;
							if(sscanf(value.c_str(),"%i,%i,%i",&r,&g,&b) == 3){
								PenArray[i].pen->SetColour(r,g,b);
							} else {
								wxString msg;
								msg.sprintf("%s.color is not a set of ints r,g,b", PenArray[i].nameabr.c_str());
								wxLogMessage(msg);
								gvpAborted = true;
							}
						}
						//style
						key.sprintf("%s.style", PenArray[i].nameabr.c_str());
						value = config[key];
						if (value != wxEmptyString) {
							//if there is anything with that Style in the map
							if (styleStringToIntMap.count(value)) {
								PenArray[i].pen->SetStyle(styleStringToIntMap[value]);
							} else {
								wxString msg;
								msg.sprintf("%s.style is not valid", PenArray[i].nameabr.c_str());
								//XXX todo: figure out a better way to log this
								cout << msg << "\n";
								//				wxLogMessage(msg);
								//				gvpAborted = true;
							}
						}
					}

					//now that the arcs and everything is initialized toggle the visibility of the nodes
					//we couldn't do it when we inited the nodes because the node visibility function
					//toggles the visibility of the its arcs, which didn't exist until just now
					//we keep this boolean to see if there are any invisible nodes and notify the user
					//just so they know
					bool contains_invisibile_nodes = false;
					if(!gvpAborted){
						for(unsigned int i = 0; i < nodes.size(); i++){
							// node visibility
							key.sprintf( "nodes[%d].visible", i);
							value = config[key];
							if (value != wxEmptyString) {
								long visible;
								if (value.ToLong(&visible) && (visible == true || visible == false)) {
									nodes[i]->setVisible(visible);
									//this indicates if there are any invisible nodes
									contains_invisibile_nodes |= !visible;
								} else {
									wxString msg;
									msg.sprintf( "nodes[%d].visible is not a number 0 or 1", i);
									wxLogMessage(msg);
									gvpAborted = true;
								}
							}
						}
					}
					//if there are invisible nodes give the user a warning just so they know and don't think
					//their data not being read in correctly
					if (!gvpAborted && contains_invisibile_nodes){
						wxString msg;
						msg.sprintf("Just so you know:\n"
								"\"%s\" contains invisible nodes\n" "if you want to make all nodes visible select "
								"View->\"Make All Nodes Visible\" from the menu tool bar.", gvpFile_cstr);
						wxLogMessage(msg);
					}
					// At this point we no longer need the file parser, and since it
					// uses global variables, we can't have more than one. Thus, we
					// delete it now, rather than later.
					delete fp;
					fp = NULL;

					// mangle the name and set the tab title and status bar
					int beginning, ending;
					if ( (beginning = file.rfind('/') + 1) < 1 )
						beginning = 0;
					ending = file.rfind(wxT(".str"));
					if (ending <= beginning)
						ending = file.rfind(wxT(".gvp"));
					if (ending <= beginning)
						ending = file.size();
					//actually make the page since we parsed the data
					parentNotebook->InsertPage(parentNotebook->GetPageCount(),
							this, file.substr(beginning, ending-beginning),
							true);
					//indicate that we parsed sucessfully
					data_parsed = true;
				} else {
					//exited on error or crashed
					wxString msg;
					msg.sprintf(
							"Could not parse the strFile: \"%s\"\n"
							"so we cannot load it, check the terminal output for more data",
							strFile_cstr);
					wxLogMessage(msg);
				}
			}
		} else {
			//we don't have a valid strfile
			wxString msg;
			msg.sprintf(
					"Position file:\n"
					"\"%s\"\n"
					"claims that its strfile is:"
					"\n \"%s\"\n"
					"but that file does not exist or cannot be opened.\n"
					"Please edit the position file to point to a valid\n"
					"strfile and re-open the Position File.", gvpFile_cstr, strFile_cstr);
			wxLogMessage(msg);
		}
	}
	fp = NULL;

	if (data_parsed){
		// Update the status bar and stuff
		onComeForward();
		// Item 7 in gZoomMap is a 1.0 scaling factor.
		setScale(ZOOM_1_INDEX);
	}

	//we haven't done a save as yet so there is no old gvp file
	gvpFile_old = "";
}

/**
 *******************************************************************
 * Deletes all the dynamically allocated data including frame
 * separators, arcs, nodes, and the drawing buffer.
 *
 * \pre Should not be called explicitly.
 *
 * \post The StructPage object will have freed all its dynamically
 *	  allocated memory and become invalid.
 *
 * \note The StructPage object will have freed all its dynamically
 *	  allocated memory and become invalid.
 *
 * \return Nothing.
 *******************************************************************/
StructPage::~StructPage( void )
{
	int numNodes = nodes.size();

	// delete all the frameEnds
	for (int i = frameEnds.size() - 1; i >= 0; i--) {
		delete frameEnds[i];
		frameEnds[i] = NULL;
	}
	// delete all the frameNameTags
	for (int i = frameNameTags.size() - 1; i >= 0; i--) {
		delete frameNameTags[i];
		frameNameTags[i] = NULL;
	}
	// delete all the arcs
	for (int i = numNodes - 1; i >= 0; i--) {
		for (int j = numNodes - 1; j >= 0; j--) {
			if (arcs[i][j]) {
				delete arcs[i][j];
				arcs[i][j] = NULL;
			}
		}
	}
	// delete all the nodes
	for (int i = numNodes - 1; i >= 0; i--) {
		delete nodes[i];
		nodes[i] = NULL;
	}
	// delete the drawing buffer
	if (content != NULL)
		delete content;
	content = NULL;
}

/**
 *******************************************************************
 * Fills the config map with the keys and values specified in the
 * position file.
 *
 * \pre gvpFile must reference a valid position file.
 *
 * \post Any keys and values specified in the position file will be
 *	  put in the config map.
 *
 * \note Any keys and values specified in the position file will be
 *	  put in the config map.
 *
 * \return void
 *******************************************************************/
void
StructPage::fillMap( void )
{
	wxTextFile gvp(gvpFile);
	wxString line;

	if (gvp.Open()) {
		// get the keys and values
		if (gvp.GetLineCount()) {
			for ( line = gvp.GetFirstLine();
					!gvp.Eof();
					line = gvp.GetNextLine() ) {
				config[line.BeforeFirst('=')] = line.AfterFirst('=');
			}
			config[line.BeforeFirst('=')] = line.AfterFirst('=');
		}
		// get the one absolutely required key=value pair
		strFile = config[wxT("strFile")];
		// complain if it wasn't present
		if (!strFile.length()) {
			wxLogMessage(wxT("gvp file doesn't specify a structure file"));
			gvpAborted = true;
		} 
	} else {
		wxLogMessage(wxString::Format("Failed to open position file \"%s\" for reading", gvpFile.c_str()));
		gvpAborted = true;
	}
}

/**
 *******************************************************************
 * Places the nodes it finds in the structure file in their grid
 * positions or in the position specified in the config map (from
 * the position file) if one was specified.
 *
 * \pre The FileParser should be valid and have already parsed the
 *	  structure file.
 *
 * \post The vector of nodes will be populated with positioned nodes.
 *
 * \note The vector of nodes will be populated with positioned nodes.
 *
 * \return void
 *******************************************************************/
void
StructPage::initNodes( void )
{
	unsigned int curFrame = 0;
	wxPoint curPos;

	curPos.x = 200*ACTUAL_SCALE;
	curPos.y = 120*ACTUAL_SCALE;
	nodes.clear();
	
	int numVars = fp->numVarsInPrologue + fp->numVarsInChunk +
	fp->numVarsInEpilogue;
	firstChunkFrame = fp->_firstChunkframe;
	firstEpilogueFrame = fp->_lastChunkframe + 1;
	numFrames = fp->_maxFrame + 1;
	int numRows = (int)ceil(sqrt(2.0*numVars));
	int yMax = 0;
	wxString key, value;

	std::map< RVInfo::rvParent, unsigned int > nameGvpNodeMap;

	// check that gvp has valid info and matches str
	key.sprintf("numNodes");
	value = config[key];
	if (value != wxEmptyString) {
		long numNodes;
		if (value.ToLong(&numNodes)) {
			if (numNodes != numVars) {
				wxLogMessage("gvp and str disagree on number of nodes");
				gvpAborted = true;
			}
			// map rvId's to the number in the gvp file
			RVInfo::rvParent rvId;
			for (int i = 0; i < numNodes; i++) {
				key.sprintf("nodes[%d].nametag.name", i);
				value = config[key];
				if (value != wxEmptyString) {
					rvId.first = value;
					key.sprintf("nodes[%d].frame", i);
					value = config[key];
					if (value != wxEmptyString) {
						long temp;
						if (value.ToLong(&temp)) {
							rvId.second = temp;
							nameGvpNodeMap[rvId] = i;
						} else {
							wxString msg;
							msg.sprintf("nodes[%d].frame is not a number", i);
							wxLogMessage(msg);
							gvpAborted = true;
						}
					}
				}
			}
		} else {
			wxLogMessage("numNodes is not a number");
			gvpAborted = true;
		}
	}
	// check that gvp has valid info and matches str
	key.sprintf("numFrames");
	value = config[key];
	if (value != wxEmptyString) {
	long gvpNumFrames;
	if (value.ToLong(&gvpNumFrames)) {
		if (gvpNumFrames != (long)numFrames) {
			wxLogMessage("gvp and str disagree on number of frames");
			gvpAborted = true;
		}
	} else {
		wxLogMessage("numFrames is not a number");
		gvpAborted = true;
	}
	}

	// for each node
	firstNodeInFrame.push_back(nodes.size()); // should be 0
	assert(firstNodeInFrame.size() == curFrame+1);
	for (int i = 0, row = 0; i < numVars; i++, row++) {
		/* If this is a new frame or there are too many nodes in this
		 * column, then move over to the next column. */
		if (row >= numRows || fp->rvInfoVector[i].frame != curFrame) {
			if (curPos.y > yMax)
				yMax = curPos.y;
			curPos.x += 80*ACTUAL_SCALE;
			curPos.y = 120*ACTUAL_SCALE;
			row = 0;
		}
		/* If this is a new frame, add a frame separator and move over
		 * yet again. */
		if (fp->rvInfoVector[i].frame != curFrame) {
			VizSep::FrameSepType sepType;
			if (curFrame < fp->_firstChunkframe-1)
				sepType = VizSep::PROLOGUE;
			else if (curFrame == fp->_firstChunkframe-1)
				sepType = VizSep::BEGIN_CHUNK;
			else if ( fp->_firstChunkframe <= curFrame &&
					curFrame < fp->_lastChunkframe )
				sepType = VizSep::CHUNK;
			else if (curFrame == fp->_lastChunkframe)
				sepType = VizSep::END_CHUNK;
			else //if (fp->_lastChunkframe < curFrame)
				sepType = VizSep::EPILOGUE;
			frameEnds.push_back(new VizSep( curPos.x, this, sepType ));
			// if a position was specified, move the frameEnd to it
			key.sprintf("frameEnds[%d].x", i);
			value = config[key];
			if (value != wxEmptyString) {
				long xPos;
				if (value.ToLong(&xPos)) {
					frameEnds[i]->x = xPos;
				} else {
					wxString msg;
					msg.sprintf("frameEnds[%d].x is not a number", i);
					wxLogMessage(msg);
					gvpAborted = true;
				}
			}
			// move to the next column
			curPos.x += 80*ACTUAL_SCALE;
			curFrame++;
			assert(fp->rvInfoVector[i].frame == curFrame);
			firstNodeInFrame.push_back(nodes.size());
			assert(firstNodeInFrame.size() == curFrame+1);
		}
		// move down to the next row
		curPos.y += 80*ACTUAL_SCALE;
		// add the node with this position and move it if necessary later
		nodes.push_back(new VizNode( curPos, &fp->rvInfoVector[i], this ));
		assert(nodes.size() == (unsigned)i + 1);
		nodeNameTags.push_back(&nodes[i]->nametag);
		assert(nodeNameTags.size() == nodes.size());
		// add it to the map from Id to index
		nameVizNodeMap[nodes[i]->rvId] = i;

		// if a position was specified for a variable with the same
		// name, move the node to it
		// x
		key.sprintf( "nodes[%d].center.x",
				 nameGvpNodeMap.count(nodes[i]->rvId) ?
				 nameGvpNodeMap[nodes[i]->rvId] : i );
		value = config[key];
		if (value != wxEmptyString) {
			long xPos;
			if (value.ToLong(&xPos)) {
				nodes[i]->center.x = xPos;
				nodeNameTags[i]->pos.x = xPos;
			} else {
				wxString msg;
				msg.sprintf( "nodes[%d].center.x is not a number",
					nameGvpNodeMap.count(nodes[i]->rvId) ?
					nameGvpNodeMap[nodes[i]->rvId] : i );
				wxLogMessage(msg);
				gvpAborted = true;
			}
		}
		// y
		key.sprintf( "nodes[%d].center.y",
				 nameGvpNodeMap.count(nodes[i]->rvId) ?
				 nameGvpNodeMap[nodes[i]->rvId] : i );
		value = config[key];
		if (value != wxEmptyString) {
			long yPos;
			if (value.ToLong(&yPos)) {
				nodes[i]->center.y = yPos;
				nodeNameTags[i]->pos.y = yPos;
			} else {
				wxString msg;
				msg.sprintf( "nodes[%d].center.y is not a number",
					nameGvpNodeMap.count(nodes[i]->rvId) ?
					nameGvpNodeMap[nodes[i]->rvId] : i );
				wxLogMessage(msg);
				gvpAborted = true;
			}
		}

		// and move the nametag, if requested
		// x
		key.sprintf( "nodes[%d].nametag.pos.x",
				 nameGvpNodeMap.count(nodes[i]->rvId) ?
				 nameGvpNodeMap[nodes[i]->rvId] : i );
		value = config[key];
		if (value != wxEmptyString) {
			long xPos;
			if (value.ToLong(&xPos)) {
				nodeNameTags[i]->pos.x = xPos;
			} else {
				wxString msg;
				msg.sprintf( "nodes[%d].nametag.pos.x is not a number",
					nameGvpNodeMap.count(nodes[i]->rvId) ?
					nameGvpNodeMap[nodes[i]->rvId] : i );
				wxLogMessage(msg);
				gvpAborted = true;
			}
		}
		// y
		key.sprintf( "nodes[%d].nametag.pos.y",
				 nameGvpNodeMap.count(nodes[i]->rvId) ?
				 nameGvpNodeMap[nodes[i]->rvId] : i );
		value = config[key];
		if (value != wxEmptyString) {
			long yPos;
			if (value.ToLong(&yPos)) {
				nodeNameTags[i]->pos.y = yPos;
			} else {
				wxString msg;
				msg.sprintf( "nodes[%d].nametag.pos.y is not a number",
					nameGvpNodeMap.count(nodes[i]->rvId) ?
					nameGvpNodeMap[nodes[i]->rvId] : i );
				wxLogMessage(msg);
				gvpAborted = true;
			}
		}
		// nametag visibility
		key.sprintf( "nodes[%d].nametag.visible",
				 nameGvpNodeMap.count(nodes[i]->rvId) ?
				 nameGvpNodeMap[nodes[i]->rvId] : i );
		value = config[key];
		if (value != wxEmptyString) {
			long visible;
			if (value.ToLong(&visible) && (visible == true || visible == false)) {
				nodeNameTags[i]->setVisible(visible);
			} else {
				wxString msg;
				msg.sprintf( "nodes[%d].nametag.visible is not a number 0 or 1",
					nameGvpNodeMap.count(nodes[i]->rvId) ?
					nameGvpNodeMap[nodes[i]->rvId] : i );
				wxLogMessage(msg);
				gvpAborted = true;
			}
		}
	}
	firstNodeInFrame.push_back(nodes.size());
	assert(firstNodeInFrame.size() == (unsigned)numFrames+1);
	// get an idea of how big the canvas should be
	if (curPos.y > yMax)
		yMax = curPos.y;
	canvasWidth = curPos.x + 200*ACTUAL_SCALE;
	canvasHeight = yMax + 200*ACTUAL_SCALE;
	// use the width from the gvp if one was specified
	key.sprintf("canvasWidth");
	value = config[key];
	if (value != wxEmptyString) {
		long w;
		if (value.ToLong(&w)) {
			canvasWidth = w;
		} else {
			wxString msg;
			msg.sprintf("canvasWidth is not a number");
			wxLogMessage(msg);
			gvpAborted = true;
		}
	}
	// use the height from the gvp if one was specified
	key.sprintf("canvasHeight");
	value = config[key];
	if (value != wxEmptyString) {
		long h;
		if (value.ToLong(&h)) {
			canvasHeight = h;
		} else {
			wxString msg;
			msg.sprintf("canvasHeight is not a number");
			wxLogMessage(msg);
			gvpAborted = true;
		}
	}
	
	
	/* Since I'm not sure what we made up and what was in the file,
	 * consider the file dirty. */
	gvpDirty = true;
	// Update the status bar to reflect that we're dirty.
	onComeForward();
}

/**
 *******************************************************************
 * Creates the arcs and control points defined in the position
 * file plus the mandatory beginning and ending points, or makes
 * them up if they weren't specified.
 *
 * \pre The nodes should already be created and placed.
 *
 * \post The arcs and control points are created for all the arcs in
 *	  the structure file according to the position file.
 *
 * \note The arcs and control points are created for all the arcs in
 *	  the structure file according to the position file.
 *
 * \return void
 *******************************************************************/
void
StructPage::initArcs( void )
{
	unsigned int numNodes = fp->rvInfoVector.size();
	std::vector< VizArc* > curVec;

	// initialize an empty arc matrix
	arcs.clear();
	for (unsigned int i = 0; i < numNodes; i++) {
		curVec.clear();
		arcs.push_back(curVec);
		assert(arcs.size() == i + 1);
		for (unsigned int j = 0; j < numNodes; j++) {
			arcs[i].push_back(NULL);
			assert(arcs[i].size() == j + 1);
		}
	}

	// for each node
	for (unsigned int j = 0; j < numNodes; j++) {
		// for each switching parent
		for ( unsigned int k = 0;
			  k < nodes[j]->rvi->switchingParents.size();
			  k++ ) {
			/* The Id in the parent vector is relative to the current
			 * frame, so I need to translate it to an absolute frame
			 * number given a zero unrolling. */
			RVInfo::rvParent absId;
			absId.first = nodes[j]->rvi->switchingParents[k].first;
			absId.second = nodes[j]->rvi->switchingParents[k].second
			+ nodes[j]->rvi->frame;
			// figure out which node number that is
			int i = nameVizNodeMap[ absId ];
			// if the arc doesn't exist yet, make it
			if (!arcs[i][j])
				arcs[i][j] = newArc(i, j);
			// label it as a switching arc (an arc to a switching parent)
			arcs[i][j]->setSwitching();
			
		}
		// for each conditional parent
		for ( unsigned int k = 0;
				k < nodes[j]->rvi->conditionalParents.size();
				k++ ) {
			for ( unsigned int l = 0;
					l < nodes[j]->rvi->conditionalParents[k].size();
					l++ ) {
				/* The Id in the parent vector is relative to the
				 * current frame, so I need to translate it to an
				 * absolute frame number given a zero unrolling. */
				RVInfo::rvParent absId;
				absId.first = nodes[j]->rvi->conditionalParents[k][l].first;
				absId.second = nodes[j]->rvi->conditionalParents[k][l].second
					+ nodes[j]->rvi->frame;
				// figure out which node number that is
				int i=nameVizNodeMap[absId];
				// if the arc doesn't exist yet, make it
				if (!arcs[i][j])
					arcs[i][j] = newArc(i, j);
				/* label it as random or deterministic */
				if (nodes[j]->rvi->rvType == RVInfo::t_discrete){
						if(nodes[j]->rvi->discImplementations[k] == CPT::di_MTCPT)
							arcs[i][j]->setDeterministic();
						else
							arcs[i][j]->setRand();
				} else
					arcs[i][j]->setRand();
			}
		}
	}
	/** \todo Only mark the file dirty when it really is, maybe. On
	 *  the other hand, the file could have just about anything in it
	 *  to start with, and those things won't necessarily get written
	 *  back. Then again, that may be a reason to not mark the file
	 *  dirty, since saving won't affect the graph but will only get
	 *  rid of whatever other info was in the file. */
	gvpDirty = true;
	// Update the status bar to reflect that we're dirty.
	onComeForward();
}

/**
 *******************************************************************
 * Given two node indexes, create and return an arc from the
 * first node to the second. Take the definition from the config
 * map (from the position file) if possible. Otherwise make up a
 * couple control points for it.
 *
 * \param i The index of the first ("from") node in the nodes vector.
 * \param j The index of the second ("to") node in the nodes vector.
 *
 * \pre The nodes need to exist and be initialized.
 *
 * \post A VizArc will be dynamically allocated and returned to you
 *	  along with appropriate ControlPoints.
 *
 * \note This method dynamically allocates memory that the caller is
 *	  expected to dispose of.
 *
 * \return A VizArc will be dynamically allocated and returned to you
 *	  along with appropriate ControlPoints.
 *******************************************************************/
VizArc*
StructPage::newArc( int i, int j )
{
	std::vector< ControlPoint* > *cps
	= new std::vector< ControlPoint* >;
	wxString key, value;

	cps->clear();
	// The first control point _must_ be the first node's center point.
	cps->push_back(new ControlPoint(nodes[i]->center));

	wxPoint pt;
	/* Check to see if the arc is in the config map and if so how many
	 * control points it should have. */
	key.sprintf("arcs[%d][%d].numCPs", i, j);
	value = config[key];
	if (value != wxEmptyString) {
		long numCPs;
		if (value.ToLong(&numCPs)) {
			// for each pointrol point in the arc
			for (int k = 1; k < numCPs - 1; k++) {
				// try to get the x position
				key.sprintf("arcs[%d][%d].cps[%d].pos.x", i, j, k);
				value = config[key];
				if (value != wxEmptyString) {
					long xPos;
					if (value.ToLong(&xPos)) {
						// and try to get the y position
						key.sprintf("arcs[%d][%d].cps[%d].pos.y", i, j, k);
						value = config[key];
						if (value != wxEmptyString) {
							long yPos;
							/* If you get good values for x and y,
							 * then add the control point. */
							if (value.ToLong(&yPos)) {
								pt.x = xPos;
								pt.y = yPos;
								cps->push_back(new ControlPoint(pt));
							} else {
								wxString msg;
								msg.sprintf("arcs[%d][%d].cps[%d].pos.y "
										"is not a number", i, j, k);
								wxLogMessage(msg);
								gvpAborted = true;
							}
						}
					} else {
						wxString msg;
						msg.sprintf("arcs[%d][%d].cps[%d].pos.x "
								"is not a number", i, j, k);
						wxLogMessage(msg);
						gvpAborted = true;
					}
				}
			}
		} else {
			wxString msg;
			msg.sprintf("arcs[%d][%d].numCPs is not a number", i, j);
			wxLogMessage(msg);
			gvpAborted = true;
		}
	} else {
		/* Make up a couple points so that arcs don't overlap and hide
		 * each other. */
		pt.x = (2*nodes[i]->center.x + nodes[j]->center.x)/3 + 20*ACTUAL_SCALE;
		pt.y = (2*nodes[i]->center.y + nodes[j]->center.y)/3 - 20*ACTUAL_SCALE;
		ControlPoint *mid = new ControlPoint(pt);
		cps->push_back(mid);
		pt.x = (nodes[i]->center.x + 2*nodes[j]->center.x)/3 + 20*ACTUAL_SCALE;
		pt.y = (nodes[i]->center.y + 2*nodes[j]->center.y)/3 + 20*ACTUAL_SCALE;
		mid = new ControlPoint(pt);
		cps->push_back(mid);
	}

	/* The last control point absolutely must be the center of the
	 * second node. */
	cps->push_back(new ControlPoint(nodes[j]->center));
	return new VizArc(cps, this);
}

/**
 *******************************************************************
 * Event handler for key presses that simply deletes control points
 * when the delete key is pressed.
 *
 * \param event The event generated by the keypress. At the very
 *	  least, event.m_keyCode should be filled in.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Any selected ControlPoints will be deleted iff the event says
 *	  that the delete key was pressed.
 *
 * \note Any selected ControlPoints will be deleted iff the event says
 *	  that the delete key was pressed.
 *
 * \return void
 *******************************************************************/
void
StructPage::OnChar( wxKeyEvent &event )
{
	//get the mouse position incase we need it
	wxPoint mouse_pos = wxGetMousePosition();
	ScreenToClient(&mouse_pos.x,&mouse_pos.y);
	CalcUnscrolledPosition( mouse_pos.x, mouse_pos.y, &mouse_pos.x, &mouse_pos.y );
	// calculate unscaled position
	mouse_pos.x = (int)round(mouse_pos.x / gZoomMap[displayScale]);
	mouse_pos.y = (int)round(mouse_pos.y / gZoomMap[displayScale]);

	if (event.m_keyCode == WXK_DELETE || event.m_keyCode == 'r') {
		save_undo();
		//if we actually deleted any cps then redraw, otherwise pop
		//one state from the undo stack
		if (deleteSelectedCps()){
			redraw();
			blit();
		} else
			restore_state(undo_stack);
	} else if (event.m_keyCode == WXK_LEFT) {
		save_undo();
		moveSelected(-ACTUAL_SCALE*(event.ShiftDown()?10:1), 0);
		redraw();
		blit();
	} else if (event.m_keyCode == WXK_UP) {
		save_undo();
		moveSelected(0, -ACTUAL_SCALE*(event.ShiftDown()?10:1));
		redraw();
		blit();
	} else if (event.m_keyCode == WXK_RIGHT) {
		save_undo();
		moveSelected(ACTUAL_SCALE*(event.ShiftDown()?10:1), 0);
		redraw();
		blit();
	} else if (event.m_keyCode == WXK_DOWN) {
		save_undo();
		moveSelected(0, ACTUAL_SCALE*(event.ShiftDown()?10:1));
		redraw();
		blit();
	} else if (event.m_keyCode == 's'){
		//toggle nametag's visiblity
		//find a node under the mouse
		int node_index = nodeUnderPt(mouse_pos);
		//if we find a node then toggle it's nametag visibility
		if (0 <= node_index){
			if (nodes[node_index]->isVisible()){
				save_undo();
				nodes[node_index]->nametag.toggleVisible();
				redraw();
				blit();
			}
		} else {
			//see if we're on a frame name tag
			int tag_index = -1;
			for (unsigned int i = 0; i < frameNameTags.size(); i++) {
				if (frameNameTags[i]->onMe(mouse_pos)){
					tag_index = i;
					break;
				}
			}
			//if this is true then we're on a frameNameTag
			if (0 <= tag_index){
				frameNameTags[tag_index]->toggleVisible();
				redraw();
				blit();
			}
		}
	} else if (event.m_keyCode == 'w'){
		bool something_toggled = false;
		save_undo();
		//toggle the visiblity of the selected nodes and framenametags
		for (int i = 0; i < (int)nodes.size(); i++) {
			if (nodes[i]->getSelected() && nodes[i]->isVisible()){
				nodes[i]->nametag.toggleVisible();
				something_toggled = true;
			}
		}
		for (int i = 0; i < (int)frameNameTags.size(); i++) {
			if (frameNameTags[i]->getSelected()){
				frameNameTags[i]->toggleVisible();
				something_toggled = true;
			}
		}
		if(something_toggled){
			redraw();
			blit();
		} else
			restore_state(undo_stack);
	} else if (event.m_keyCode == 'a'){
		popUpNodeInfo(mouse_pos);
	} else if (event.m_keyCode == 'd'){
		//toggle node and its edges visiblity
		//find a node under the mouse
		int node_index = nodeUnderPt(mouse_pos);
		//if there is a node under the mouse, toggle it
		if(node_index >= 0 && node_index < (int)nodes.size()){
			save_undo();
			toggleNodeAndArcVisibility(node_index);
			redraw();
			blit();
		}
	} else if (event.m_keyCode == 'e'){
		save_undo();
		//toggle visibility for all the nodes in the current selection
		bool node_toggled = false;
		for (int i = 0; i < (int)nodes.size(); i++) {
			if (nodes[i]->getSelected()){
				node_toggled = true;
				toggleNodeAndArcVisibility(i);
			}
		}
		//if no node was toggled then we don't need to redraw, and we should pop
		//from the undo stack
		if (node_toggled){
			redraw();
			blit();
		} else
			restore_state(undo_stack);
	} else if (event.m_keyCode == 'f'){
		//highlight
		//find a node under the mouse
		int node_index = nodeUnderPt(mouse_pos);
		nextHighlightState(node_index);

		redraw();
		blit();
	} else if (event.GetKeyCode() == 'z' || event.GetKeyCode() == 'Z'){
		//select current frame
		//if it is capitol then add it to the current selection
		int frame_n = 0;
		if(event.GetKeyCode() != 'Z')
			setAllSelected(false);
		for (int i = frameEnds.size() - 1; i >= 0; i--) {
			if(mouse_pos.x > frameEnds[i]->x){
				frame_n = i + 1;
				break;
			}
		}
		selectAllInFrame(frame_n);
		redraw();
		blit();
	} else if(event.GetKeyCode() == 'x'){
		//invert selection
		for (unsigned int i = 0; i < nodes.size(); i++){
			nodes[i]->setSelected(!nodes[i]->getSelected());
		}
		redraw();
		blit();
	} else if(event.GetKeyCode() == 'g'){
		if (undo()){
			setScale(displayScale);
			redraw();
			blit();
			MainVizWindow->UpdateMenuChecks(this);
		} else {
			//make beep sound
			printf("\a");
			fflush(stdout);
		}
	} else if(event.GetKeyCode() == 'G'){
		if (redo()){
			setScale(displayScale);
			redraw();
			blit();
			MainVizWindow->UpdateMenuChecks(this);
		} else {
			//make beep sound
			printf("\a");
			fflush(stdout);
		}
	} else {
		event.Skip();
	}
}

/**
 *******************************************************************
 * Given a node index (which may be -1 to indicate that we're not on
 * a node), toggle the Node and its arcs visibility
 *
 * \pre The StructPage must be fully initialized as well as its nodes and arcs
 *
 * \post the node and its incoming and outgoing arcs will have their visibility
 * toggled, if the node (node_index) was highlighted everything will be unhighlighted
 *
 * \note the node and its incoming and outgoing arcs will have their visibility
 * toggled, if the node (node_index) was highlighted everything will be unhighlighted
 *
 * \return void
 *******************************************************************/

void
StructPage::toggleNodeAndArcVisibility(int node_index){
	//first make sure it is in range
	int numNodes = nodes.size();
	if (0 <= node_index && node_index < (int)nodes.size()){
		nodes[node_index]->toggleVisible();
		//if this node is highlighed make all it's parents and edges
		//unhighlighted
		if (nodes[node_index]->isHighlighted()){
			for (int i = 0; i < numNodes; i++) {
				for (int j = 0; j < numNodes; j++) {
					if (arcs[i][j])
						arcs[i][j]->setHighlighted(false);
					nodes[i]->setHighlightState(off);
				}
			}
			//unhighlight the node
			nodes[node_index]->setHighlightState(off);
		}
	}
}

/**
 *******************************************************************
 * Given a node index (which may be -1 to indicate that we're not on
 * a node), go to the next Highlighting state.  If we're not on a node
 * then unhighlight everything.
 *
 * \pre The StructPage must be fully initialized.
 *
 * \post The highlighting state for the given node will be changed,
 * it's parents, children and in/out arcs may also be highlighted or
 * unhighlighted depending on the state of the given node.  If a negative
 * node id is given then all nodes and arcs will be set to an unhighlighted
 * state
 *
 * \note The highlighting state for the given node will be changed,
 * it's parents, children and in/out arcs may also be highlighted or
 * unhighlighted depending on the state of the given node.  If a negative
 * node id is given then all nodes and arcs will be set to an unhighlighted
 * state
 *
 * \return void
 *******************************************************************/

void
StructPage::nextHighlightState(int node_index){
	//highlight
	int numNodes = nodes.size();
	highlight_states next_highlight_state = off;
	//if it is a valid node get it's highlight state
	if (0 <= node_index && node_index < numNodes){
		//if the node isn't visible then we shouldn't highlight anything
		if (!nodes[node_index]->isVisible())
			node_index = -1;
		else{
			//figure out the next state for highlighting
			if(off == nodes[node_index]->getHighlightState())
				next_highlight_state = children;
			else
				next_highlight_state = nodes[node_index]->getNextHighlightState();
		}
	}

	switch (next_highlight_state){
		case children:
			//step through all the arcs, if it starts at our node
			//highlight it, also highlight the destination node
			for(int i = 0; i < numNodes; i++){
				//turn off highlighting for all but our node's children
				if (arcs[node_index][i])
					nodes[i]->setHighlightState(on);
				else
					nodes[i]->setHighlightState(off);

				if ( i == node_index){
					for (int j = 0; j < numNodes; j++) {
						//j is a child of our node
						if (arcs[node_index][j])
							arcs[node_index][j]->setHighlighted(true);
					}
				} else {
					for (int j = 0; j < numNodes; j++) {
						if (arcs[i][j])
							arcs[i][j]->setHighlighted(false);
					}
				}
			}
			nodes[node_index]->setHighlightState(children);
			break;
		case parents:
			//step though all the arcs, if it ends at our node
			//highlight it, also highlight the source node (parent)
			for(int i = 0; i < numNodes; i++){
				//turn off all highlighting for nodes, though we may turn it
				//back on later
				nodes[i]->setHighlightState(off);
				for (int j = 0; j < numNodes; j++) {
					if (arcs[i][j]) {
						if (j == node_index){
							//i is a parent
							nodes[i]->setHighlightState(on);
							arcs[i][j]->setHighlighted(true);
						} else 
							arcs[i][j]->setHighlighted(false);
					} 
				}
			}
			nodes[node_index]->setHighlightState(parents);
			break;
		case both:
			//step through all the arcs, if it starts or ends at our node
			//highlight it and the nodes at either end
			for(int i = 0; i < numNodes; i++){
				for (int j = 0; j < numNodes; j++) {
					if (arcs[i][j]) {
						if (i == node_index){
							//i is a parent
							nodes[j]->setHighlightState(on);
							arcs[i][j]->setHighlighted(true);
						} else if (j == node_index){
							//i is a child
							nodes[i]->setHighlightState(on);
							arcs[i][j]->setHighlighted(true);
						} else {
							arcs[i][j]->setHighlighted(false);
						}
					} 
				}
			}
			nodes[node_index]->setHighlightState(both);
			break;
		default:
			//clear everything
			for(int i = 0; i < numNodes; i++){
				nodes[i]->setHighlightState(off);
				for (int j = 0; j < numNodes; j++) {
					if (arcs[i][j])
						arcs[i][j]->setHighlighted(false);
				}
			}
			break;
	}
}

/**
 *******************************************************************
 * Iterate through all of the Nodes, if they are in the given frame
 * and visible then select them, also select their control points
 *
 * \pre The StructPage must be fully initialized.
 *
 * \post All the visible nodes in the given frame will be selected.
 *
 * \note All the visible nodes in the given frame will be selected.
 *
 * \return void
 *******************************************************************/

void
StructPage::selectAllInFrame(int frame)
{
	int numNodes = nodes.size();
	//make sure the frame number is valid
	if (frame < 0 || frame >= numFrames)
		return;
	//for each node, if it is in the frame and is visible select it,
	//also select the controlpoints along its arcs
	for (int i = 0; i < numNodes; i++){
		if (nodes[i]->rvId.second == frame && nodes[i]->isVisible()){
			nodes[i]->setSelected(true);
			for (int j = 0; j < numNodes; j++){
				if (arcs[i][j]){
					for(unsigned int k = 1; k < arcs[i][j]->cps->size() - 1; k++){
						(*arcs[i][j]->cps)[k]->setSelected(true);
					}
				}
				if (arcs[j][i]){
					for(unsigned int k = 1; k < arcs[j][i]->cps->size() - 1; k++){
						(*arcs[j][i]->cps)[k]->setSelected(true);
					}
				}
			}
		}
	}
	
}

/**
 *******************************************************************
 * Iterate through all of the ControlPoints (except those that
 * correspond to the beginning or end of an arc and are therefore at
 * the center of a node) and delete any which are currently selected.
 *
 * \pre The StructPage must be fully initialized.
 *
 * \post Any selected ControlPoints will be deleted.
 *
 * \note Any selected ControlPoints will be deleted.
 *
 * \return bool, if any control points have been deleted
 *******************************************************************/
bool
StructPage::deleteSelectedCps( void )
{
	int numNodes = nodes.size();
	bool cps_deleted = false;
	// iterate through every possible arc
	for (int i = 0; i < numNodes; i++) {
		for (int j = 0; j < numNodes; j++) {
			if (arcs[i][j]) {
				int end = arcs[i][j]->cps->size() - 1;
				// and every control point
				for (int k = end - 1; k > 0; k--) {
					// delete only if it was selected
					if ((*arcs[i][j]->cps)[k]->getSelected()) {
						// delete's on it's own?
						arcs[i][j]->cps->erase(arcs[i][j]->cps->begin() + k);
						arcs[i][j]->points->DeleteNode(arcs[i][j]->points->Item(k));
						cps_deleted = true;
					}
				}
			}
		}
	}
	/* Assume that this made it dirty. We could make this accurate at
	 * the cost of an assignment per deletion. */
	gvpDirty = true;
	// Update the status bar to indicate that we're dirty.
	onComeForward();
	return cps_deleted;
}

/**
 *******************************************************************
 * Handles mouse events. This is where all the action takes
 * place. This keeps track of the selection and moves things around.
 *
 * \param event Describes the details of the mouse event, such as
 * where it occurred, what modifier keys are held down, and what kind
 * of event it is.
 *
 * \pre event should be properly filled out and the StructPage should
 * be fully initialized.
 *
 * \post Any number of things might have changed since this is the
 * main driver for user manipulation of the graph.
 *
 * \note Numerous side effects.
 *
 * \return void
 *******************************************************************/
void
StructPage::OnMouseEvent( wxMouseEvent &event )
{
	static bool dragging = false;
	static bool boxSelecting = false;
	static bool shifted = false;
	static bool start_drag = false;
	static wxPoint dragStart;
	bool screenDirty = false;
	wxPoint pt;

	event.GetPosition(&pt.x, &pt.y);
	CalcUnscrolledPosition( pt.x, pt.y, &pt.x, &pt.y );
	// calculate unscaled position
	pt.x = (int)round(pt.x / gZoomMap[displayScale]);
	pt.y = (int)round(pt.y / gZoomMap[displayScale]);
	
	// find out what (if anything) was clicked
	Selectable *pointee = itemAt(pt);

	if (event.LeftDown()) {
		shifted = event.ShiftDown(); // keep this around
		if (pointee) {
			if (shifted) {
				// LeftDown + Shifted + on something
				pointee->toggleSelected();
				screenDirty = true;
			} else {
				// LeftDown + not Shifted + on something
				dragging = true;
				start_drag = true;	//we are starting a drag
				dragStart.x = pt.x;
				dragStart.y = pt.y;
				if (!pointee->getSelected()) {
					setAllSelected(false);
					pointee->setSelected(true);
					screenDirty = true;
				}
			}
		} else {
			// LeftDown + not on anything
			boxSelecting = true;
			selectBox.x = pt.x;
			selectBox.y = pt.y;
			selectBox.width = 0;
			selectBox.height = 0;
			drawSelectBox = true;
			if (!shifted) {
				// LeftDown + not Shifted + not on anything
				setAllSelected( false );
				screenDirty = true;
			}
		}
	} else if (event.Dragging() && event.LeftIsDown()) {
		if (boxSelecting) {
			// LeftDragging + box selecting
			selectBox.width = pt.x - selectBox.x;
			selectBox.height = pt.y - selectBox.y;
		} else if (dragging) {
			//if this is the beg of a drag then we should save the undo
			if(start_drag){
				save_undo();
				start_drag = false;
			}
			// LeftDragging + moving items
			moveSelected( pt.x - dragStart.x, pt.y - dragStart.y );
			dragStart.x = pt.x;
			dragStart.y = pt.y;
		}
		screenDirty = true;
	} else if (event.LeftUp()) {
		start_drag = false;
		if (boxSelecting) {
			// LeftUp + box selecting = done box selecting
			if (pt.x < selectBox.x) {
				// swap pt.x and selectBox.x so we have a real rectangle
				int temp = pt.x;
				pt.x = selectBox.x;
				selectBox.x = temp;
			}
			if (pt.y < selectBox.y) {
				// swap pt.y and selectBox.y so we have a real rectangle
				int temp = pt.y;
				pt.y = selectBox.y;
				selectBox.y = temp;
			}
			selectBox.width = pt.x - selectBox.x;
			selectBox.height = pt.y - selectBox.y;
			if (!shifted) {
				// LeftUp + box selecting + not shifted when we started
				setAllSelected(false);
			}
			toggleSelectedInRect(selectBox);
			screenDirty = true;
		} else if (dragging && snapToGrid) {
			// LeftUp + dragging + snap to grid = do grid snapping
			snapSelectedToGrid();
			screenDirty = true;
		}
		// Whatever we were doing, we're done.
		boxSelecting = shifted = dragging = drawSelectBox = false;
	} else if (event.RightDown()) {
		// clone a control point maybe
		ControlPoint *cp = dynamic_cast<ControlPoint *>(pointee);
		/* If the user right clicked on a control point that makes our
		 * job a bit easier. */
		if (cp) {
			save_undo();
			int index = -1;
			// get the arc and control point index
			VizArc *arc = findArcOwning(cp, index);
			if (arc && index >= 0) {
				// make a new one
				wxPoint where(pt.x - NEW_CP_OFFSET, pt.y - NEW_CP_OFFSET);
				arc->cps->insert( arc->cps->begin() + index,
						new ControlPoint(where) );
				// Tell it who its parent is.
				(*arc->cps)[index]->arc = arc;
				// must keep these in sync
				arc->points->Insert( index,
						(wxObject*)&(*arc->cps)[index]->pos );
				screenDirty = true;
			}
		} else {
			// check to see if it's on a line between two nodes
			// that have no other control points
			int numNodes = nodes.size();
			VizArc *arc = NULL;
			for (int i = 0; i < numNodes && !arc; i++) {
				for (int j = 0; j < numNodes && !arc; j++) {
					arc = arcs[i][j];
					if (arc == NULL) continue;
					/* This stuff is all used to figure out if a point
					 * is near the line. */
					int x0 = (*arc->cps)[0]->pos.x, x1 = (*arc->cps)[1]->pos.x,
						 y0 = (*arc->cps)[0]->pos.y, y1 = (*arc->cps)[1]->pos.y;
					double y = ( (x1-x0) ?
							(y1-y0)/(double)(x1-x0)*(pt.x-x0) + y0 :
							(y0+y1)/2 ),
							 x = ( (y1-y0) ?
									 (x1-x0)/(double)(y1-y0)*(pt.y-y0) + x0 :
									 (x0+x1)/2 ),
							 epsilon = fabs(2*ACTUAL_SCALE*(y1-y0)/(double)(x1-x0)),
							 delta = fabs(2*ACTUAL_SCALE*(x1-x0)/(double)(y1-y0));
					epsilon = (epsilon>2*ACTUAL_SCALE?epsilon:2*ACTUAL_SCALE);
					delta = ( delta>2*ACTUAL_SCALE ? delta : 2*ACTUAL_SCALE );
#if 0
					// allow this to compile for some insight
					wxString msg;
					msg.sprintf("(pt.x, pt.y) = (%d, %d)\n"
							"(x0, y0)	 = (%d, %d)\n"
							"(x1, y1)	 = (%d, %d)\n"
							"(x, y)	   = (%f, %f)\n"
							"(del, eps)   = (%f, %f)\n"
							"fabs(x-pt.x) = %f\n"
							"fabs(y-pt.y) = %f\n",
							pt.x, pt.y, x0, y0, x1, y1, x, y,
							delta, epsilon, fabs(x-pt.x), fabs(y-pt.y) );
					wxLogMessage(msg);
#endif
					// if the arc exists and has only two control points
					if ( arc->cps->size() == 2 &&
							// and the click is in the rectangle
							( ( x0-1 <= pt.x && pt.x <= x1+1 ) ||
							  ( x1-1 <= pt.x && pt.x <= x0+1 ) ) &&
							( ( y0-1 <= pt.y && pt.y <= y1+1 ) ||
							  ( y1-1 <= pt.y && pt.y <= y0+1 ) ) &&
							// and (roughly) on the line
							fabs(y-pt.y) <= epsilon && fabs(x-pt.x) <= delta ) {
						//save undo
						save_undo();
						// make a new control point
						wxPoint where(pt.x-NEW_CP_OFFSET, pt.y-NEW_CP_OFFSET);
						arc->cps->insert( arc->cps->begin() + 1,
								new ControlPoint(where) );
						// Tell it who its parent is.
						(*arc->cps)[1]->arc = arc;
						// must keep these in sync
						arc->points->Insert(1,(wxObject*)&(*arc->cps)[1]->pos);
						screenDirty = true;
					} else {
						arc = NULL;
					}
				}
			}
		}
	} else if (event.MiddleDown()){
		popUpNodeInfo(pt);
	}
	if ( screenDirty ) {
		redraw();
		blit();
	}

	// Anything else to be done?
	event.Skip();
}

/**
 *******************************************************************
 * Pops up a dialog for User Customizations of Pens
 *
 * \param Void
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post 1 or More pens may be customized, the StructPage will be redrawn
 *
 * \note Numerous side effects.
 *
 * \return void
 *******************************************************************/
void
StructPage::OnCustomizePen()
{
	save_undo();
	GmtkPenSelectDialog * pen_dialog = new GmtkPenSelectDialog(this, -1,
			wxT("Pen Customization Dialog"), PenArray, numPens);
	if(pen_dialog->ShowModal() != wxID_OK){
		undo();
	} else {
		redraw();
		blit();
	}
}

/**
 *******************************************************************
 * If the point is on a node then it pops up information about that Variable
 *
 * \param a Point
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post a temporary window will pop up that can be closed by clicking anywhere
 *
 * \note Numerous side effects.
 *
 * \return void
 *******************************************************************/

void
StructPage::popUpNodeInfo(wxPoint pt){
	//print information about a node and it's parents
	wxString message;
	int i = nodeUnderPt(pt);	//find a node under point pt
	//if we're on a node
	if(0 <= i) {
		//name
		message = "Name: ";
		message.append(nodes[i]->nametag.name);
		message.append("\n");
		message.append(wxString::Format("Frame: %d\n", nodes[i]->rvId.second));
		//type
		message.append("Type: ");
		switch (nodes[i]->rvi->rvType){
			case RVInfo::t_discrete:
				message.append("discrete\n");
				break;
			case RVInfo::t_continuous:
				message.append("continuous\n");
				break;
			default:
				message.append("unknown\n");
		}

		//disposition
		message.append("Disposition: ");
		switch (nodes[i]->rvi->rvType){
			case RVInfo::d_hidden:
				message.append("hidden\n");
				break;
			case RVInfo::d_observed:
				message.append("observed\n");
				break;
			default:
				message.append("unknown\n");
		}

		//cardinality
		message.append(wxString::Format("Cardinality: %d\n",nodes[i]->rvi->rvCard));

		//switching parents
		if (nodes[i]->rvi->switchingParents.size() > 0){
			for (unsigned int j = 0; j < nodes[i]->rvi->switchingParents.size(); j++){
				message.append("Switching Parent: ");
				message.append(nodes[i]->rvi->switchingParents[j].first.c_str());
				message.append("(");
				message.append(wxString::Format("%d", nodes[i]->rvi->switchingParents[j].second));
				message.append(")");
				message.append(" using ");

				//mappings
				message.append(" mapping(\"");
				switch(nodes[i]->rvi->switchMapping.liType){
					case RVInfo::ListIndex::li_String :
						message.append(nodes[i]->rvi->switchMapping.nameIndex.c_str());
						break;
					case RVInfo::ListIndex::li_Index :
						message.append(wxString::Format("%d", nodes[i]->rvi->switchMapping.intIndex));
						break;
					default:
						message.append("unknown");
				}
				message.append("\")\n");
			}
		} else
			message.append("Switching Parents: nil\n");

		//conditional parents
		if (nodes[i]->rvi->conditionalParents.size() > 0){
			message.append("Conditional Parents: ");
			for (unsigned int j = 0; j < nodes[i]->rvi->conditionalParents.size(); j++){
				if (nodes[i]->rvi->conditionalParents[j].size() == 0)
					message.append("nil");
				for (unsigned int k = 0; k < nodes[i]->rvi->conditionalParents[j].size(); k++){
					message.append(nodes[i]->rvi->conditionalParents[j][k].first.c_str());
					message.append("(");
					message.append(wxString::Format("%d", nodes[i]->rvi->conditionalParents[j][k].second));
					message.append(")");
					if (k != nodes[i]->rvi->conditionalParents[j].size() - 1)
						message.append(", ");
				}
				message.append(" using ");
				//discrete implementations
				if (nodes[i]->rvi->rvType == RVInfo::t_discrete){
					switch (nodes[i]->rvi->discImplementations[j]){
						case CPT::di_MDCPT: 
							message.append("DenseCPT");
							break;
						case CPT::di_MSCPT: 
							message.append("SparseCPT");
							break;
						case CPT::di_MTCPT: 
							message.append("DeterministicCPT");
							break;
						case CPT::di_USCPT: 
							message.append("UnityScoreCPT");
							break;
						case CPT::di_NGramCPT: 
							message.append("Ngram\"languagemodel\"CPT");
							break;
						case CPT::di_FNGramCPT: 
							message.append("factoredngram\"languagemodel\"CPT");
							break;
						case CPT::di_VECPT:
							message.append("VirtualEvidenceCPT");
							break;
						case CPT::di_LatticeNodeCPT:
							message.append("latticenodeCPT");
							break;
						case CPT::di_LatticeEdgeCPT:	
							message.append("latticeedgeCPT");
							break;
						default:
							message.append("unknown");
					}
					//continuous implementations
				} else if (nodes[i]->rvi->rvType == RVInfo::t_continuous){
					switch (nodes[i]->rvi->contImplementations[j]){ 
						case MixtureCommon::ci_unknown:
							message.append("unknown");
							break;
						case MixtureCommon::ci_mixture:
							message.append("mixture ");
							break;
						case MixtureCommon::ci_gausSwitchMixture:
							message.append("gausSwitchMixGaussian");
							break;
						case MixtureCommon::ci_logitSwitchMixture:
							message.append("logitSwitchMixGaussian");
							break;
						case MixtureCommon::ci_mlpSwitchMixture:
							message.append("mlpSwitchMixGaussian");
							break;
						case MixtureCommon::ci_zeroScoreMixture:
							message.append("zeroScoreMixGaussian");
							break;
						case MixtureCommon::ci_unityScoreMixture:
							message.append("unityScoreMixGaussian");
							break;
						default:
							message.append("unknown");
					};
				} else {
					message.append("unknown");
				}
				//listIndices
				if(nodes[i]->rvi->rvType == RVInfo::t_continuous){
					if(nodes[i]->rvi->listIndices[j].collectionName != ""){
						message.append("collection(\"");
						message.append(nodes[i]->rvi->listIndices[j].collectionName.c_str());
						message.append("\") ");
					}
					message.append("mapping");
				}
				switch(nodes[i]->rvi->listIndices[j].liType){
					case RVInfo::ListIndex::li_String :
						message.append("(\"");
						message.append(nodes[i]->rvi->listIndices[j].nameIndex.c_str());
						message.append("\")");
						break;
					case RVInfo::ListIndex::li_Index :
						message.append(wxString::Format("(\"%d\")", nodes[i]->rvi->listIndices[j].intIndex));
						break;
					default:
						message.append("(unknown)");
				}

				//if there are more parents [based on switching varaibles], print | between them
				if (j != nodes[i]->rvi->conditionalParents.size() - 1)
					message.append(" | ");
			}
		} else
			message.append("Conditional Parents: nil");
		new wxTipWindow(this, message, 500, NULL);
	}
}


/**
 *******************************************************************
 * Determine the parent arc and index in that arc for the given
 * ControlPoint.
 *
 * \param cp The ControlPoint whose parent and index you want.
 * \param index A reference to an integer where the index into the
 * vector of control points is to be stored.
 *
 * \pre \p cp must exist and be properly initialized.
 *
 * \post \p index holds the index into the control points vector or -1
 * if it couldn't be found.
 *
 * \note No side effects.
 *
 * \return A pointer to the VizArc that owns the given control point.
 *******************************************************************/
VizArc*
StructPage::findArcOwning( ControlPoint *cp, int& index )
{
	/* Originally, ControlPoints didn't have arc pointers, so this
	 * involved iterating through all the control points to find a
	 * match. Now we just check the arc pointer assuming that it is
	 * accurate. */
	if (cp->arc) {
		int end = cp->arc->cps->size() - 1;
		for (int k = 1; k < end; k++) {
			/* If this is the right control point, set the index and
			 * return the arc. */
			if ( (*cp->arc->cps)[k] == cp ) {
				index = k;
				return cp->arc;
			}
		}
	}
	index = -1;
	return NULL; // nothing found
}

/**
 *******************************************************************
 * Convenience method to blit the drawing buffer to the screen.
 *
 * \pre The StructPage must be fully initialized. In particular, we
 * need the drawing buffer to be valid.
 *
 * \post The screen will be updated.
 *
 * \note The screen will be updated.
 *
 * \return void
 *******************************************************************/
void
StructPage::blit( void )
{
	wxClientDC dc( this );
	PrepareDC( dc );
	blit( dc );
}

/**
 *******************************************************************
 * Copy the contents of the drawing buffer to the specified device
 * context.
 *
 * \pre The StructPage must be fully initialized. In particular, we
 * need the drawing buffer to be valid.
 *
 * \post The given device context will be updated.
 *
 * \note The given device context will be updated.
 *
 * \return void
 *******************************************************************/
void
StructPage::blit( wxDC& dc )
{
	// Clear away what used to be there to keep the background clean
	dc.SetBackground(*wxLIGHT_GREY_BRUSH);
	wxCoord w, h;
	dc.GetSize(&w, &h);
	wxRegion tempClip(0, 0, w, h);
	tempClip.Subtract(wxRegion( 0, 0,
		(int)round(getWidth()*gZoomMap[displayScale]),
		(int)round(getHeight()*gZoomMap[displayScale]) ));
	dc.SetClippingRegion(tempClip);
	dc.Clear();
	dc.DestroyClippingRegion();
	// all the action is the constructor and destructor
	wxBufferedDC bdc( &dc, *content );
}

/**
 *******************************************************************
 * Determines whether a node (or its nametag) should prevent moving a
 * frame separator from \p x0 to \p x1.
 *
 * \param x0 The x position the frame separator is moving from.
 * \param x1 The x position the frame separator is moving to.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post No real change.
 *
 * \note No side effects.
 *
 * \return \c true if the movement would cross a node or its nametag
 * and should therefor be prevented, or \c false if the move should be
 * fine.
 *******************************************************************/
bool
StructPage::crossesNode( wxCoord x0, wxCoord x1 )
{
	int numNodes = nodes.size();

	for (int i = 0; i < numNodes; i++) {
		if ( ( ( x0 < nodes[i]->center.x-NODE_RADIUS !=
				x1 < nodes[i]->center.x-NODE_RADIUS ) ||
					( x0 > nodes[i]->center.x+NODE_RADIUS !=
				x1 > nodes[i]->center.x+NODE_RADIUS ) ) ||
					( ( x0 < nodes[i]->nametag.pos.x !=
				x1 < nodes[i]->nametag.pos.x ) ||
					( x0 > nodes[i]->nametag.pos.x + 2*NODE_RADIUS !=
				x1 > nodes[i]->nametag.pos.x + 2*NODE_RADIUS ) ) )
			return true;
	}
	for (unsigned int i = 0; i < frameNameTags.size(); i++) {
		if ( ( x0 < frameNameTags[i]->pos.x !=
				x1 < frameNameTags[i]->pos.x ) ||
			 ( x0 > frameNameTags[i]->pos.x + 2*NODE_RADIUS !=
				x1 > frameNameTags[i]->pos.x + 2*NODE_RADIUS ) )
			return true;
	}
	return false;
}


/**
 *******************************************************************
 * Determines when moving an item (node or nametag) from x position
 * \p x_cur to x position \p x_new move out of its frame
 *
 * \pre The StructPage (and in particular the vector of frame
 * separators) should be fully initialized.
 *
 * \post No real changes.
 *
 * \note No side effects.
 *
 * \return \c true if moving in bounds of the frame or moving toward being in bounds
 *	\c false if moving out of bounds
 *******************************************************************/
bool
StructPage::itemMovingInBounds(int frame, wxCoord x_cur, wxCoord x_new)
{
	int nframe_ends = frameEnds.size();
	wxCoord x_min, x_max;

	if (frame < 0){
		cout << "item has an invalid frame number\n";
		return false;
	}

	if (nframe_ends == 0){
		//there is only 1 frame
		x_min = 0;
		x_max = getWidth();
	} else if (frame == 0){
		//first frame
		x_min = 0;
		x_max = frameEnds[0]->x;
	} else if (frame == nframe_ends){
		//the item is in the last frame
		x_min = frameEnds[nframe_ends - 1]->x;
		x_max = getWidth();
	} else {
		//this item is in a frame that isn't the first or last
		x_min = frameEnds[frame - 1]->x;
		x_max = frameEnds[frame]->x;
	}

	if (x_new > x_min && x_new < x_max)
		return true;
	else {
		//is it atleast moving toward being in bounds?
		if (x_new < x_min && ((x_min - x_new) < (x_min - x_cur)))
			return true;
		else if (x_new > x_max && ((x_new - x_max) < (x_cur - x_max)))
			return true;
		else
			return false;
	}
}

/**
 *******************************************************************
 * Determines whether the given point is in the drawing area or
 * not. (If not a move to this position shouldn't be allowed.)
 *
 * \pre The StructPage should be fully initialized. In particular,
 * canvasWidth and canvasHeight should be set properly.
 *
 * \post No real changes.
 *
 * \note No side effects.
 *
 * \return \c true if the coordinates are within the canvas area, \c
 * false if they are not.
 *******************************************************************/
bool
StructPage::inBounds( wxCoord x, wxCoord y )
{
	return ( x > 0 && x < getWidth() &&
		 y > 0 && y < getHeight() );
}

/**
 *******************************************************************
 * Gives the node index for a node that is under point pt.  If there
 * is no node there then returns -1;
 * 
 * \pre The StructPage should be fully initialized.  In particular the nodes
 * should be allocated.
 *
 * \post No real changes.
 *
 * \note No side effects.
 *
 * \return node index if pt is on a node, -1 if pt is not on a node
 *
 *******************************************************************/
int
StructPage::nodeUnderPt( wxPoint pt)
{
	int numNodes = nodes.size();
	//go through the nodes
	for (int i = 0; i < numNodes; i++) {
		//if we're on a node
		if( nodes[i]->onMe(pt))
			return i;
	}
	return -1;	//if we're not on a node
}

/**
 *******************************************************************
 * Move frame separator \p i \p dx units to the right, except in cases
 * where the move should not be allowed.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Frame separator \p i may be moved.
 *
 * \note Frame separator \p i may be moved.
 *
 * \return void
 *******************************************************************/
void
StructPage::moveFrameSep( int i, int dx )
{
	if ( inBounds(frameEnds[i]->x + dx, 1) &&
			!crossesNode(frameEnds[i]->x, frameEnds[i]->x + dx) ) {
		frameEnds[i]->x += dx;
	}
}

/**
 *******************************************************************
 * Move frame nametag \p i \p dx units to the right and \p dy units
 * down, except in cases where the move should not be allowed.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Frame nametag \p i may be moved.
 *
 * \note Frame nametag \p i may be moved.
 *
 * \return void
 *******************************************************************/
void
StructPage::moveFrameNameTag( int i, int dx, int dy )
{
	if ( inBounds( frameNameTags[i]->pos.x + dx,
				frameNameTags[i]->pos.y + dy ) ) {
		if(itemMovingInBounds(i, frameNameTags[i]->pos.x + NODE_RADIUS,
					frameNameTags[i]->pos.x + NODE_RADIUS + dx))
			frameNameTags[i]->pos.x += dx;
		frameNameTags[i]->pos.y += dy;
	}
}

/**
 *******************************************************************
 * Move node nametag \p i \p dx units to the right and \p dy units
 * down, except in cases where the move should not be allowed.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Node nametag \p i may be moved.
 *
 * \note Node nametag \p i may be moved.
 *
 * \return void
 *******************************************************************/
void
StructPage::moveNodeNameTag( int i, int dx, int dy )
{
	if ( inBounds( nodeNameTags[i]->pos.x + dx,
				nodeNameTags[i]->pos.y + dy ) ) {
		//Would be nice if the point was centered and we didn't use NODE_RADIUS
		//but since fonts can change the size of the nametag maybe this is the
		//best we can do for now
		if (itemMovingInBounds(nodeNameTags[i]->GetFrame(),
					nodeNameTags[i]->pos.x + NODE_RADIUS, nodeNameTags[i]->pos.x +
					NODE_RADIUS + dx))
			nodeNameTags[i]->pos.x += dx;
		nodeNameTags[i]->pos.y += dy;
	}
}

/**
 *******************************************************************
 * Move node \p i \p dx units to the right and \p dy units
 * down, except in cases where the move should not be allowed. This
 * also moves control points that are attached to the node.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Node \p i and attached control points may be moved.
 *
 * \note Node \p i and attached control points may be moved.
 *
 * \return void
 *******************************************************************/
void
StructPage::moveNode( int i, int dx, int dy )
{
	if ( inBounds(nodes[i]->center.x + dx, nodes[i]->center.y + dy) ) {
		if (!itemMovingInBounds(nodes[i]->GetFrame(), nodes[i]->center.x,
					nodes[i]->center.x + dx))
			dx = 0;
		nodes[i]->center.x += dx;
		nodes[i]->center.y += dy;
		/*nodes[i]->tipWin->Move( nodes[i]->center.x - NODE_RADIUS,
		  nodes[i]->center.y - NODE_RADIUS );*/
		int totalNodes = nodes.size();
		for ( int j = 0; j < totalNodes; j++ ) {
			// outgoing arcs
			if (arcs[i][j]) {
				(*arcs[i][j]->cps)[0]->pos.x += dx;
				(*arcs[i][j]->cps)[0]->pos.y += dy;
			}
			// incoming arcs
			if (arcs[j][i]) {
				(*arcs[j][i]->cps)[ arcs[j][i]->cps->size()-1 ]->pos.x += dx;
				(*arcs[j][i]->cps)[ arcs[j][i]->cps->size()-1 ]->pos.y += dy;
			}
		}
	}
}

/**
 *******************************************************************
 * Move control point \p k on the arc from node \p i to node \p j \p
 * dx units to the right and \p dy units down, except in cases where
 * the move should not be allowed.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Control point \p i, \p j, \p k may be moved.
 *
 * \note Control point \p i, \p j, \p k may be moved. This does not
 *   pay attention to which control point it is moving. In particular,
 *   it doesn't know if the control point should remain attached to a
 *   node.
 *
 * \return void
 *******************************************************************/
void
StructPage::moveControlPoint( int i, int j, int k, int dx, int dy )
{
	if ( inBounds( (*arcs[i][j]->cps)[k]->pos.x+dx,
		   (*arcs[i][j]->cps)[k]->pos.y+dy ) ) {
	(*arcs[i][j]->cps)[k]->pos.x += dx;
	(*arcs[i][j]->cps)[k]->pos.y += dy;
	}
}

/**
 *******************************************************************
 * Move all selected items \p dx units to the right and \p dy units
 * down, except in cases where the move should not be allowed.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Selected items may be moved.
 *
 * \note Selected items may be moved.
 *
 * \return void
 *******************************************************************/
void
StructPage::moveSelected( int dx, int dy )
{
	int numNodes = nodes.size();

	// frame separators
	for (unsigned int i = 0; i < frameEnds.size(); i++) {
		if ( frameEnds[i]->getSelected() )
			moveFrameSep( i, dx );
	}
	// frame name tags
	for (unsigned int i = 0; i < frameNameTags.size(); i++) {
		if ( frameNameTags[i]->getSelected() && frameNameTags[i]->isVisible())
			moveFrameNameTag(i, dx, dy);
	}
	// node name tags
	for (int i = 0; i < numNodes; i++) {
		if ( nodeNameTags[i]->getSelected() && nodeNameTags[i]->isVisible())
			moveNodeNameTag(i, dx, dy);
	}
	// nodes
	for (int i = 0; i < numNodes; i++) {
		//don't move the node unless it is visible
		if ( nodes[i]->getSelected() && nodes[i]->isVisible())
			moveNode(i, dx, dy);
	}
	// then arcs
	for (int i = 0; i < numNodes; i++) {
		for (int j = 0; j < numNodes; j++) {
			if (arcs[i][j]) {
				int end = arcs[i][j]->cps->size() - 1;
				assert(end >= 1);
				for ( int k = 1; k < end; k++ ) {
					if ( (*arcs[i][j]->cps)[k]->getSelected() && (*arcs[i][j]->cps)[k]->isVisible())
						moveControlPoint(i, j, k, dx, dy);
				}
			}
		}
	}
	// Things are almost definitely dirty.
	gvpDirty = true;
	// Update the status bar to reflect the dirt.
	onComeForward();
}

/**
 *******************************************************************
 * Snap all selected items to the grids.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Selected items may be moved.
 *
 * \note Selected items may be moved.
 *
 * \return void
 *******************************************************************/
void
StructPage::snapSelectedToGrid( void )
{
	int numNodes = nodes.size();

	// frame separators
	int dx, dy;
	for (unsigned int i = 0; i < frameEnds.size(); i++) {
		dx = -(frameEnds[i]->x % GRID_SIZE);
		if (dx <= -GRID_SIZE/2)
			dx += GRID_SIZE;
		if ( frameEnds[i]->getSelected() )
			moveFrameSep(i, dx);
	}
	// frame name tags
	for (unsigned int i = 0; i < frameNameTags.size(); i++) {
		dx = -(frameNameTags[i]->pos.x % GRID_SIZE);
		if (dx <= -GRID_SIZE/2)
			dx += GRID_SIZE;
		dy = -(frameNameTags[i]->pos.y % GRID_SIZE);
		if (dy <= -GRID_SIZE/2)
			dy += GRID_SIZE;
		if ( frameNameTags[i]->getSelected() )
			moveFrameNameTag(i, dx, dy);
	}
	// node name tags
	for (int i = 0; i < numNodes; i++) {
		dx = -(nodeNameTags[i]->pos.x % GRID_SIZE);
		if (dx <= -GRID_SIZE/2)
			dx += GRID_SIZE;
		dy = -(nodeNameTags[i]->pos.y % GRID_SIZE);
		if (dy <= -GRID_SIZE/2)
			dy += GRID_SIZE;
		if ( nodeNameTags[i]->getSelected() )
			moveNodeNameTag(i, dx, dy);
	}
	// nodes
	for (int i = 0; i < numNodes; i++) {
		dx = -(nodes[i]->center.x % GRID_SIZE);
		if (dx <= -GRID_SIZE/2)
			dx += GRID_SIZE;
		dy = -(nodes[i]->center.y % GRID_SIZE);
		if (dy <= -GRID_SIZE/2)
			dy += GRID_SIZE;
		if ( nodes[i]->getSelected() )
			moveNode(i, dx, dy);
	}
	// then arcs
	for (int i = 0; i < numNodes; i++) {
		for (int j = 0; j < numNodes; j++) {
			if (arcs[i][j]) {
				int end = arcs[i][j]->cps->size();
				assert(end > 1);
				for ( int k = 0; k < end; k++ ) {
					dx = -((*arcs[i][j]->cps)[k]->pos.x % GRID_SIZE);
					if (dx <= -GRID_SIZE/2)
						dx += GRID_SIZE;
					dy = -((*arcs[i][j]->cps)[k]->pos.y % GRID_SIZE);
					if (dy <= -GRID_SIZE/2)
						dy += GRID_SIZE;
					if ( (*arcs[i][j]->cps)[k]->getSelected() )
						moveControlPoint(i, j, k, dx, dy);
				}
			}
		}
	}
	// Things are almost definitely dirty.
	gvpDirty = true;
	// Update the status bar to reflect the dirt.
	onComeForward();
}

/**
 *******************************************************************
 * Copy the layout from frame \p from to frame \p to (at least as far
 * as they have common variables).
 *
 * \param from The frame from which to copy the node layout.
 * \param to The frame to which to copy the node layout.
 *
 * \pre The StructPage should be fully initialized and \p from and \to
 *   should be valid frame numbers.
 *
 * \post The nodes in frame \p to may have been moved to positions
 *   similar to those in the \p from frame.
 *
 * \note The nodes in frame \p to may have been moved to positions
 *   similar to those in the \p from frame.
 *
 * \return void
 *******************************************************************/
void
StructPage::copyFrameLayout( int from, int to )
{
	RVInfo::rvParent destNodeId;
	destNodeId.second = to;
	int destNode, dx, dxFromLeft, dxFromRight, dy, iLeft, iRight, destLeft,
		 destRight, iOffsetFromLeft, iOffsetFromRight, destOffsetFromLeft,
		 destOffsetFromRight;
	for (int i = firstNodeInFrame[from]; i < firstNodeInFrame[from+1]; i++) {
		destNodeId.first = nodes[i]->rvId.first;
		if (nameVizNodeMap.count(destNodeId)) {
			destNode = nameVizNodeMap[destNodeId];

			// Figure out where to move the node.
			iLeft = ( from==0 ? 0 : frameEnds[from-1]->x );
			iOffsetFromLeft = nodes[i]->center.x - iLeft;
			iRight = ( from==numFrames-1 ? getWidth() : frameEnds[from]->x );
			iOffsetFromRight = nodes[i]->center.x - iRight;

			destLeft = ( to==0 ? 0 : frameEnds[to-1]->x );
			destOffsetFromLeft = nodes[destNode]->center.x - destLeft;
			destRight = to==numFrames-1 ? getWidth() : frameEnds[to]->x;
			destOffsetFromRight = nodes[destNode]->center.x - destRight;

			dxFromLeft = iOffsetFromLeft - destOffsetFromLeft;
			dxFromRight = iOffsetFromRight - destOffsetFromRight;

			dx = abs(dxFromLeft)<=abs(dxFromRight) ? dxFromLeft : dxFromRight;
			dy = nodes[i]->center.y - nodes[destNode]->center.y;
			moveNode(destNode, dx, dy);

			// Figure out where to move the node's nametag.
			iLeft = ( from==0 ? 0 : frameEnds[from-1]->x );
			iOffsetFromLeft = nodeNameTags[i]->pos.x - iLeft;
			iRight = ( from==numFrames-1 ? getWidth() : frameEnds[from]->x );
			iOffsetFromRight = nodeNameTags[i]->pos.x - iRight;

			destLeft = ( to==0 ? 0 : frameEnds[to-1]->x );
			destOffsetFromLeft = nodeNameTags[destNode]->pos.x - destLeft;
			destRight = to==numFrames-1 ? getWidth() : frameEnds[to]->x;
			destOffsetFromRight = nodeNameTags[destNode]->pos.x - destRight;

			dxFromLeft = iOffsetFromLeft - destOffsetFromLeft;
			dxFromRight = iOffsetFromRight - destOffsetFromRight;

			dx = abs(dxFromLeft)<=abs(dxFromRight) ? dxFromLeft : dxFromRight;
			dy = nodeNameTags[i]->pos.y - nodeNameTags[destNode]->pos.y;
			moveNodeNameTag(destNode, dx, dy);

			// Rebuild the VizArc with similar ControlPoints
			VizArc *iArc = NULL, *destArc = NULL;
			RVInfo::rvParent destParentId, destChildId;
			int destParent, destChild, numNodes = nodes.size();
			for (int j = 0; j < numNodes; j++) {
				iArc = arcs[i][j];
				// if the arc doesn't exist, we're done for now
				if (!iArc)
					continue;
				destChildId.first = nodes[j]->rvId.first;
				destChildId.second = nodes[j]->rvId.second -
					nodes[i]->rvId.second + nodes[destNode]->rvId.second;
				// if there is no similar node in relation to
				// destNode, then we're done for now
				if (!nameVizNodeMap.count(destChildId))
					continue;
				destChild = nameVizNodeMap[destChildId];
				destArc = arcs[destNode][destChild];
				// if that arc doesn't exist we're done for now
				if (!destArc)
					continue;
				// We've found two similar arcs. Now we need to
				// make destArc look like iArc.
				copyArcLayout(i, j, destNode, destChild, false);
			}
			for (int j = 0; j < numNodes; j++) {
				iArc = arcs[j][i];
				// if the arc doesn't exist, we're done for now
				if (!iArc)
					continue;
				destParentId.first = nodes[j]->rvId.first;
				destParentId.second = nodes[j]->rvId.second -
					nodes[i]->rvId.second + nodes[destNode]->rvId.second;
				// if there is no similar node in relation to
				// destNode, then we're done for now
				if (!nameVizNodeMap.count(destParentId))
					continue;
				destParent = nameVizNodeMap[destParentId];
				destArc = arcs[destParent][destNode];
				// if that arc doesn't exist we're done for now
				if (!destArc)
					continue;
				// We've found two similar arcs. Now we need to
				// make destArc look like iArc.
				copyArcLayout(j, i, destParent, destNode, true);
			}
		}
	}

	// also move the frame's nametag
	iLeft = ( from==0 ? 0 : frameEnds[from-1]->x );
	iOffsetFromLeft = frameNameTags[from]->pos.x - iLeft;
	iRight = ( from==numFrames-1 ? getWidth() : frameEnds[from]->x );
	iOffsetFromRight = frameNameTags[from]->pos.x - iRight;

	destLeft = ( to==0 ? 0 : frameEnds[to-1]->x );
	destOffsetFromLeft = frameNameTags[to]->pos.x - destLeft;
	destRight = to==numFrames-1 ? getWidth() : frameEnds[to]->x;
	destOffsetFromRight = frameNameTags[to]->pos.x - destRight;

	dxFromLeft = iOffsetFromLeft - destOffsetFromLeft;
	dxFromRight = iOffsetFromRight - destOffsetFromRight;

	dx = abs(dxFromLeft)<=abs(dxFromRight) ? dxFromLeft : dxFromRight;
	dy = frameNameTags[from]->pos.y - frameNameTags[to]->pos.y;
	moveFrameNameTag(to, dx, dy);

	redraw();
	blit();
}

/* this class represents a dialog that shows 2 lists and has an ok and cancel button
 * the first list allows for 1 selection and the second list allows for many selections
 * this allows for a 1 to many relationship to be specified (like copying the layout of
 * frame to many other frames
 */
class wxgmtk1toManyDialog : public wxDialog {
	public:
		wxgmtk1toManyDialog( wxWindow *parent,
				wxWindowID id,
				const wxString &title,
				wxString choices_1[],
				wxString choices_Many[],
				int choices_1_size,
				int choices_Many_size,
				wxString caption_1,
				wxString caption_Many,
				const wxPoint& position = wxDefaultPosition,
				const wxSize& size = wxDefaultSize,
				long style = wxDEFAULT_DIALOG_STYLE );
		~wxgmtk1toManyDialog();
		int Get_1_Selection() { return listbox1->GetSelection();}
		int Get_Many_Selection(wxArrayInt *selected) { return listboxMany->GetSelections(*selected);}
	private:
		void OnOk( wxCommandEvent &event );
		void OnCancel( wxCommandEvent &event );
		void OnApply( wxCommandEvent &event );
		wxListBox * listbox1;
		wxListBox * listboxMany;
		void (*function)(wxgmtk1toManyDialog);
		DECLARE_EVENT_TABLE()

};

BEGIN_EVENT_TABLE(wxgmtk1toManyDialog,wxDialog)
    EVT_BUTTON( wxID_OK, wxgmtk1toManyDialog::OnOk )
    EVT_BUTTON( wxID_CANCEL, wxgmtk1toManyDialog::OnCancel )
    EVT_BUTTON( wxID_APPLY, wxgmtk1toManyDialog::OnApply )
END_EVENT_TABLE()

wxgmtk1toManyDialog::wxgmtk1toManyDialog( wxWindow *parent,
		wxWindowID id,
		const wxString &title,
		wxString choices_1[],
		wxString choices_Many[],
		int choices_1_size,
		int choices_Many_size,
		wxString caption_1,
		wxString caption_Many,
		const wxPoint& position,
		const wxSize& size,
		long style) : 
	wxDialog( parent, id, title, position, size, style )
{
	wxBoxSizer *topsizer = new wxBoxSizer( wxVERTICAL );
	wxBoxSizer *selection_sizer = new wxBoxSizer( wxHORIZONTAL );

	wxStaticText * caption_1_ob = new wxStaticText(this,-1, caption_1, wxDefaultPosition);
	selection_sizer->Add( caption_1_ob, 0, wxALL | wxALIGN_CENTER,5 );

	listbox1 = new wxListBox((wxWindow*)this, -1, wxDefaultPosition,
			wxDefaultSize, choices_1_size, choices_1, wxLB_SINGLE | wxLB_NEEDED_SB,
			wxDefaultValidator, caption_1);
	selection_sizer->Add( listbox1, 1, wxALIGN_LEFT | wxEXPAND | wxALL, 5);
	listbox1->SetSelection(0);	//select the first item by default
	
	wxStaticText * caption_Many_op = new wxStaticText(this,-1, caption_Many, wxDefaultPosition);
	selection_sizer->Add( caption_Many_op, 0, wxALL | wxALIGN_CENTER ,5);

	listboxMany = new wxListBox((wxWindow*)this, -1, wxDefaultPosition,
			wxDefaultSize, choices_Many_size, choices_Many, wxLB_EXTENDED | wxLB_NEEDED_SB,
			wxDefaultValidator, caption_Many);
	selection_sizer->Add( listboxMany, 1, wxALIGN_LEFT | wxEXPAND | wxALL, 5);

	//add the selections, align them along the center of the dialog box
	topsizer->Add( selection_sizer, 1, wxALIGN_CENTER | wxEXPAND);

	//create buttons and add them
	wxBoxSizer *button_sizer = new wxBoxSizer( wxHORIZONTAL );
	button_sizer->Add( new wxButton( this, wxID_OK, "Ok" ), 0, wxALL, 10 );
//	button_sizer->Add( new wxButton( this, wxID_APPLY, "Apply" ), 0, wxALL, 10 );
	button_sizer->Add( new wxButton( this, wxID_CANCEL, "Cancel" ), 0, wxALL, 10 );
	//add the buttons, align them along the center of the dialog box
	topsizer->Add( button_sizer, 0, wxALIGN_CENTER );

	SetAutoLayout( true );     // tell dialog to use sizer
	SetSizer( topsizer );      // actually set the sizer

	topsizer->Fit( this );            // set size to minimum size as calculated by the sizer
	topsizer->SetSizeHints( this );   // set size hints to honour mininum size
}

wxgmtk1toManyDialog::~wxgmtk1toManyDialog(){
	delete listboxMany;
	delete listbox1;
}

void wxgmtk1toManyDialog::OnOk( wxCommandEvent &event ){
	event.Skip();
}
void wxgmtk1toManyDialog::OnCancel( wxCommandEvent &event ){
	listbox1->Clear();
	listboxMany->Clear();
	event.Skip();
}
void wxgmtk1toManyDialog::OnApply( wxCommandEvent &event ){
	event.Skip();
}


/**
 *******************************************************************
 * Ask the user which frame's layout they would like copied to which
 * and pass the buck.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Nodes may be moved around.
 *
 * \note Nodes may be moved around.
 *
 * \return void
 *******************************************************************/
void
StructPage::copyFrameLayout( void )
{
	int from = 0, count;
	wxString temp_string;
	wxArrayInt selections_Many;
	wxString *choices = new wxString [numFrames];
	int to;

	for (int i = 0; i < numFrames; i++){
		temp_string.Printf("%i", i);
		choices[i] = temp_string;
	}

	wxgmtk1toManyDialog *copyDialog = new wxgmtk1toManyDialog((wxWindow*)this,
			-1,
			wxString("Copy Frame Layout"),
			choices, choices, numFrames, numFrames,
			wxString("Copy Layout From:"),
			wxString("To:")
			);
	copyDialog->ShowModal();
	from = copyDialog->Get_1_Selection();
	
	//make sure a source has been selected and that there is at least one destination
	if(from < 0 || 1 > (count = copyDialog->Get_Many_Selection(&selections_Many))){
		delete [] choices;
		delete copyDialog;
		return;	//zero selected or canceled
	}

	save_undo();
	
	//copy the frame layout to all the selected frames (except the source if it is selected)
	for(int i = 0; i < count; i++){
		to = selections_Many.Item(i);
		if (to != from)
			copyFrameLayout(from, to);
	}
	delete [] choices;
	delete copyDialog;
}

/**
 *******************************************************************
 * Rebuild the arc from \p iTo to \p jTo based on the ControlPoints in
 * the arc from \p iFrom to jFrom (using relative
 * offsets from ends which we don't move).
 *
 * \param iFrom The beginning node of the arc to copy the layout from.
 * \param jFrom The ending node of the arc to copy the layout from.
 * \param iTo The beginning node of the arc to copy the layout to.
 * \param jTo The ending node of the arc to copy the layout to.
 * \param backward Whether to copy using relative offsets from the
 * child (rather than from the parent).
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post ControlPoints may be added, deleted, or moved.
 *
 * \note ControlPoints may be added, deleted, or moved.
 *
 * \return void
 *******************************************************************/
void
StructPage::copyArcLayout( int iFrom, int jFrom,
			   int iTo, int jTo,
			   bool backward )
{
	// save the first and last ControlPoints
	VizArc *from = arcs[iFrom][jFrom];
	VizArc *to = arcs[iTo][jTo];
	ControlPoint *toFirst = (*to->cps)[0];
	ControlPoint *toLast = (*to->cps)[to->cps->size()-1];
	ControlPoint *fromFirst = (*from->cps)[0];
	ControlPoint *fromLast = (*from->cps)[from->cps->size()-1];
	ControlPoint *newCP = NULL;
	wxPoint center(getWidth()/2, getHeight()/2);
	wxPoint offset, newPoint;
	// clear the list deleting all the ones we don't want
	for (unsigned int i = 1; i < to->cps->size() - 1; i++) {
		delete (*to->cps)[i];
		(*to->cps)[i] = NULL;
	}
	to->cps->clear();
	// reinsert the first ControlPoint
	to->cps->push_back(toFirst);
	// copy the ones from the template
	for (unsigned int k = 1; k < from->cps->size() - 1; k++) {
		if (backward) {
			offset.x = (*from->cps)[k]->pos.x - fromLast->pos.x;
			offset.y = (*from->cps)[k]->pos.y - fromLast->pos.y;
			newPoint.x = toLast->pos.x + offset.x;
			newPoint.y = toLast->pos.y + offset.y;
		} else {
			offset.x = (*from->cps)[k]->pos.x - fromFirst->pos.x;
			offset.y = (*from->cps)[k]->pos.y - fromFirst->pos.y;
			newPoint.x = toFirst->pos.x + offset.x;
			newPoint.y = toFirst->pos.y + offset.y;
		}
		// in case the newPoint is outside the canvas we create the
		// ControlPoint in the center of the canvas and then move it
		// into place
		newCP = new ControlPoint(center);
		// make it point to its parent arc
		newCP->arc = to;
		to->cps->push_back(newCP);
		moveControlPoint( iTo, jTo, k,
				newPoint.x - center.x,
				newPoint.y - center.y );
	}
	// put the end point on again
	to->cps->push_back(toLast);
	// now rebuild the wxList used for drawing along the points
	to->points->Clear();
	for (unsigned int k = 0; k < to->cps->size(); k++) {
		to->points->Append( (wxObject*)&(*to->cps)[k]->pos );
	}
}

/**
 *******************************************************************
 * Ask the user which partition's layout they would like copied to which
 * and pass the buck.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Nodes may be moved around.
 *
 * \note Nodes may be moved around.
 *
 * \return void
 *******************************************************************/
void
StructPage::copyPartitionLayout( void )
{
	enum {INVALID = -1, PROLOGUE, CHUNK, EPILOGUE};
	int from, to, fromStart, toStart, numFramesToCopy,
	firstFrameInFrom, lastFrameInFrom, numFramesInFrom,
	firstFrameInTo, lastFrameInTo, numFramesInTo;
	int count_Many_selected;	//number of selections
	bool backward;
	wxArrayInt selections_Many;
	wxString choices[] = {wxT("Prologue"), wxT("Chunk"), wxT("Epilogue")};
	
	wxgmtk1toManyDialog *copyDialog = new wxgmtk1toManyDialog((wxWindow*)this,
			-1,
			wxString("Copy Partition Layout"),
			choices, choices, 3, 3,
			wxString("Copy Layout From:"),
			wxString("To:")
			);
	copyDialog->ShowModal();
	from = copyDialog->Get_1_Selection();

	//make sure a source has been selected and that there is at least one destination
	if(from < 0 || 1 > (count_Many_selected = copyDialog->Get_Many_Selection(&selections_Many))){
		return;	//zero destinations selected or canceled
	}
	save_undo();
	//copy the frame layout to all the selected frames (except the source if it is selected)
	for(int i = 0; i < count_Many_selected; i++){
		to = selections_Many.Item(i);
		if (to != from){
			backward = to < from;
			firstFrameInFrom = ( from==PROLOGUE ? 0
					: ( from==CHUNK ? firstChunkFrame
						: firstEpilogueFrame ) );
			lastFrameInFrom = ( from==PROLOGUE ? firstChunkFrame - 1
					: ( from==CHUNK ? firstEpilogueFrame - 1
						: numFrames - 1 ) );
			numFramesInFrom = lastFrameInFrom - firstFrameInFrom + 1;
			firstFrameInTo = ( to==PROLOGUE ? 0
					: ( to==CHUNK ? firstChunkFrame
						: firstEpilogueFrame ) );
			lastFrameInTo = ( to==PROLOGUE ? firstChunkFrame - 1
					: ( to==CHUNK ? firstEpilogueFrame - 1
						: numFrames - 1 ) );
			numFramesInTo = lastFrameInTo - firstFrameInTo + 1;
			numFramesToCopy = ( numFramesInTo<numFramesInFrom
					? numFramesInTo
					: numFramesInFrom );
			if (backward) {
				fromStart = firstFrameInFrom;
				toStart = lastFrameInTo - numFramesToCopy + 1;
			} else {
				fromStart = lastFrameInFrom - numFramesToCopy + 1;
				toStart = firstFrameInTo;
			}
			for (int j = 0; j < numFramesToCopy; j++) {
				copyFrameLayout(fromStart + j, toStart + j);
			}
		}
	}

}

/**
 *******************************************************************
 * Return whatever item is at the given point (or NULL if there's
 * nothing there).
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post No real changes.
 *
 * \note No side effects.
 *
 * \return A pointer to an item at the given point, or \c NULL if
 * there is no such item.
 *******************************************************************/
Selectable*
StructPage::itemAt( const wxPoint& pt )
{
	// XXX: This may be the best way, but is it?
	int numNodes = nodes.size();

	if (drawCPs) {
	// then arcs
	for (int i = 0; i < numNodes; i++) {
		for (int j = 0; j < numNodes; j++) {
		if (arcs[i][j]) {
			int end = arcs[i][j]->cps->size() - 1;
			assert(end > 0);
			for ( int k = 1; k < end; k++ ) {
			if ((*arcs[i][j]->cps)[k]->onMe(pt))
				return (*arcs[i][j]->cps)[k];
			}
		}
		}
	}
	}
	if (drawFrameNames) {
	// frame name tags (only if they're drawn)
	for (unsigned int i = 0; i < frameNameTags.size(); i++) {
		if (frameNameTags[i]->onMe( pt ))
		return frameNameTags[i];
	}
	}
	if (drawNodeNames) {
	// node name tags (only if they're drawn)
	for (int i = 0; i < numNodes; i++) {
		if (nodeNameTags[i]->onMe( pt ))
		return nodeNameTags[i];
	}
	}

	if (drawNodes) {
	// nodes
	for (int i = 0; i < numNodes; i++) {
		if (nodes[i]->onMe( pt ))
		return nodes[i];
	}
	}
	if (drawFrameSeps) {
	// frame borders
	for (unsigned int i = 0; i < frameEnds.size(); i++) {
		if (frameEnds[i]->onMe( pt ))
		return frameEnds[i];
	}
	}
	// apparently nothing was found
	return NULL;
}

/**
 *******************************************************************
 * Handles a paint event by blitting the drawing buffer to the screen
 *
 * \param event Ignored.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The screen will be updated.
 *
 * \note The screen will be updated.
 *
 * \return void
 *******************************************************************/
void
StructPage::OnPaint( wxPaintEvent &WXUNUSED(event) )
{
	// Nothing has changed so we can skip the actual drawing and just blit 
	blit();
}

/**
 *******************************************************************
 * Convenience method to draw to the drawing buffer.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The drawing buffer is updated.
 *
 * \note The drawing buffer is updated.
 *
 * \return void
 *******************************************************************/
void
StructPage::redraw( void )
{
	wxMemoryDC dc;
	dc.SelectObject( *content );
	dc.Clear();
	// Scale everything appropriately.
	dc.SetUserScale(gZoomMap[displayScale], gZoomMap[displayScale]);
	// Do the actual drawing.
	draw(dc);
}

/**
 *******************************************************************
 * Draw the graph to the given device context.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The device context will be updated.
 *
 * \note The device context will be updated.
 *
 * \return void
 *******************************************************************/
void
StructPage::draw( wxDC& dc )
{
	dc.BeginDrawing();
	// Don't trust the DC defaults.
	dc.Clear();
	wxFont tempFont( labelFont.GetPointSize(), labelFont.GetFamily(),
			 labelFont.GetStyle(), labelFont.GetWeight() );
	dc.SetFont(tempFont);
	dc.SetPen(*wxBLACK_PEN);
	dc.SetTextBackground(*wxWHITE);
	dc.SetBackgroundMode(wxTRANSPARENT);
	dc.SetTextForeground(*wxBLACK);
	dc.SetBackground(*wxWHITE_BRUSH);
	dc.SetBrush(*wxWHITE_BRUSH);

	int numNodes = nodes.size();

	if (drawGrids) {
		// draw a grid to help align things
		wxPen oldPen = dc.GetPen();
		dc.SetPen(gridPen);
		for (int x = GRID_SIZE; x < getWidth(); x += GRID_SIZE) {
			dc.DrawLine(x, 0, x, getHeight());
		}
		for (int y = GRID_SIZE; y < getHeight(); y += GRID_SIZE) {
			dc.DrawLine(0, y, getWidth(), y);
		}
		dc.SetPen(oldPen);
	}

	if (drawFrameSeps) {
		// draw frame/chunk separators
		for (unsigned int i = 0; i < frameEnds.size(); i++) {
			frameEnds[i]->draw(&dc);
		}
	}

	if (drawFrameNames) {
		// draw frame labels
		for (unsigned int i = 0; i < frameNameTags.size(); i++) {
			frameNameTags[i]->draw(&dc);
		}
	}

	// then draw arcs
	for (int i = 0; i < numNodes; i++) {
		for (int j = 0; j < numNodes; j++) {
			if (arcs[i][j]) arcs[i][j]->draw(&dc, VizArc::DRAW_ARCS);
		}
	}

	if (drawNodes) {
		// then draw nodes
		wxPen oldPen = dc.GetPen();
		dc.SetPen(nodePen);
		for (int i = 0; i < numNodes; i++) {
			nodes[i]->draw(&dc);
		}
		dc.SetPen(oldPen);
	}

	if (drawFrameNames) {
		// draw frame names
		for (unsigned int i = 0; i < frameNameTags.size(); i++) {
			frameNameTags[i]->draw(&dc);
		}
	}

	if (drawNodeNames) {
		// draw node names
		for (int i = 0; i < numNodes; i++) {
			nodeNameTags[i]->draw(&dc);
		}
	}

	// then draw control points
	for (int i = 0; i < numNodes; i++) {
		for (int j = 0; j < numNodes; j++) {
			if (arcs[i][j])
				arcs[i][j]->draw(&dc, VizArc::DRAW_CPS);
		}
	}

	// then draw the selection box (maybe)
	if (drawSelectBox) {
		wxBrush oldBrush = dc.GetBrush();
		wxPen oldPen = dc.GetPen();
		dc.SetBrush(*wxTRANSPARENT_BRUSH);
		dc.SetPen(*wxGREY_PEN);
		dc.DrawRectangle(selectBox);
		dc.SetBrush(oldBrush);
		dc.SetPen(oldPen);
	}
	
	if (drawBoundingBox) {
		//draw a bounding box around the canvas
		wxPen oldPen = dc.GetPen();
		wxBrush oldBrush = dc.GetBrush();
		dc.SetPen(boundingBoxPen);
		dc.SetBrush(*wxTRANSPARENT_BRUSH);
		dc.DrawRectangle(0, 0, getWidth(), getHeight());
		dc.SetPen(oldPen);
		dc.SetBrush(oldBrush);
	}

	dc.EndDrawing();
}

/**
 *******************************************************************
 * Writes a bunch of key=value pairs to gvpFile so that they can be
 * read in again by fillMap().
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The file whose name is in gvpFile may be overwritten with the
 * current graph data.
 *
 * \note The file whose name is in gvpFile may be overwritten with the
 * current graph data.
 *
 * \return void
 *******************************************************************/
void
StructPage::Save( void )
{

	if (gvpFile.length()) {
		wxFFile gvp;
		if (gvp.Open(gvpFile, "w")) {
			wxString line;
			int numNodes = nodes.size();

			line.sprintf(".=This file was created by gmtkViz.\n");
			gvp.Write(line);

			/* 
			 * Here we run into a problem, if the strFile is in the same
			 * directory (or under) the gvp file, this line should be 
			 * relative the gvpFile location, but then the gvp file and
			 * strFile must stay in the same location relative to eachother
			 *
			 * If the strFile isn't in the same directory (or under) as the
			 * gvpFile, then we should store an absolute path for the strFile
			 * but then we cannot move the strFile without making the gvpFile's
			 * reference invalid.
			 */

			struct stat strFile_stat;
			//we need temp strings because basename and dirname alter their args
			wxString tempGVP = wxString::Format("%s",gvpFile.c_str());
			//tempFileName = dirname(gvpFile) ++ "/" ++ strFile
			wxString tempFileName = wxString::Format("%s/%s",
					dirname((char *)tempGVP.mb_str(wxConvUTF8)), 
					strFile.c_str());
			if(strFile.GetChar(0) == '/' && stat(strFile.c_str(), &strFile_stat) == 0){
				//this is an absolute path so we should use this
				line.sprintf("strFile=%s\n", strFile.c_str());
			} else if(stat(tempFileName.c_str(), &strFile_stat) == 0){
				// we check to see if the strFile is given relative to the gvpfile Location
				//this is ideal
				line.sprintf("strFile=%s\n", strFile.c_str());
			} else {
				bool strFile_error = false;
				//we need temp strings because basename and dirname alter their args
				wxString tempSTR = wxString::Format("%s",strFile.c_str());
				tempGVP = wxString::Format("%s",gvpFile.c_str());
				//see if the strFile is relative to the cwd
				tempFileName = wxString::Format("%s/%s",
						wxFileName::GetCwd().c_str(), 
						strFile.c_str());
				// if we're doing a save as, try to see if the strfile location
				// is relative to the old gvpFile
				wxString tempGVP_old = wxString::Format("%s",gvpFile_old.c_str());
				wxString tempFileName_old = wxString::Format("%s/%s/%s",
						wxFileName::GetCwd().c_str(),
						dirname((char *)tempGVP_old.mb_str(wxConvUTF8)), 
						strFile.c_str());
				if(tempSTR.Find(dirname((char *)tempGVP.mb_str(wxConvUTF8))) == 0){
					// if strfile is in directory or subdir of gvpFile
					tempGVP = wxString::Format("%s",gvpFile.c_str());
					//remove the leading directory data and the / so that strFile loc is
					//relative to gvpFile
					tempSTR.Remove(0,1 + strlen(dirname((char *)tempGVP.mb_str(wxConvUTF8))));
					line.sprintf("strFile=%s\n", tempSTR.c_str());
				} else if((gvpFile_old != "") && (stat(tempFileName_old.c_str(), &strFile_stat) == 0)){
					//strFile is relative to old gvpFile (this is a save as)
					line.sprintf("strFile=%s\n", tempFileName_old.c_str());
				} else if(stat(tempFileName.c_str(), &strFile_stat) == 0){
					//this is a absolute path because it is not in the gvpFileDir | gvpFileDir_old (or below)
					line.sprintf("strFile=%s\n", tempFileName.c_str());
				} else {
					strFile_error = true;
				}
				//if we cannot find the strFile in a resonable location give an error message
				if (strFile_error){
					wxString msg;
					msg.sprintf("strfile: \"%s\" is not in a location gmtkViz can find\n"
							"so it is likely that this gvpFile will not open correctly.\n"
							"If it doesn't open correctly then you can edit the gvpFile\n"
							"and point it to the correct location of the strFile.", strFile.c_str());
					wxMessageBox(msg,wxT("Saving"), wxOK | wxICON_ERROR);
					//just put the name of the strFile there, better than nothing
					line.sprintf("strFile=%s\n", strFile.c_str());
				}
			}
			gvp.Write(line);

			// canvas size
			line.sprintf( "canvasWidth=%d\n", getWidth());
			gvp.Write(line);
			line.sprintf( "canvasHeight=%d\n", getHeight());
			gvp.Write(line);

			// nodes (and node nametags)
			line.sprintf("numNodes=%d\n", numNodes);
			gvp.Write(line);
			for (int i = 0; i < numNodes; i++) {
				line.sprintf( "nodes[%d].center.x=%d\n",
						  i, nodes[i]->center.x );
				gvp.Write(line);
				line.sprintf( "nodes[%d].center.y=%d\n",
						  i, nodes[i]->center.y );
				gvp.Write(line);
				line.sprintf( "nodes[%d].visible=%d\n",
						  i, nodes[i]->isVisible() );
				gvp.Write(line);
				/* These are unnecessary but could be used to better
				 * match up variables and/or to verify that the gvp
				 * file goes with the str file. */
				line.sprintf( "nodes[%d].frame=%d\n",
						  i, nodes[i]->rvId.second );
				gvp.Write(line);
				line.sprintf( "nodes[%d].nametag.name=%s\n",
						  i, nodes[i]->nametag.name.c_str() );
				gvp.Write(line);
				line.sprintf( "nodes[%d].nametag.pos.x=%d\n",
						  i, nodes[i]->nametag.pos.x );
				gvp.Write(line);
				line.sprintf( "nodes[%d].nametag.pos.y=%d\n",
						  i, nodes[i]->nametag.pos.y );
				gvp.Write(line);
				line.sprintf( "nodes[%d].nametag.visible=%d\n",
						  i, nodes[i]->nametag.isVisible() );
				gvp.Write(line);
			}

			// arcs
			for (int i = 0; i < numNodes; i++) {
				for (int j = 0; j < numNodes; j++) {
					if (arcs[i][j]) {
					line.sprintf( "arcs[%d][%d]=1\n", i, j);
					gvp.Write(line);
					line.sprintf( "arcs[%d][%d].numCPs=%d\n",
							  i,j,arcs[i][j]->cps->size() );
					gvp.Write(line);
					for (unsigned int k=1;k<arcs[i][j]->cps->size()-1;k++){
						line.sprintf("arcs[%d][%d].cps[%d]"
							 ".pos.x=%d\n", i, j, k,
							 (*arcs[i][j]->cps)[k]->pos.x);
						gvp.Write(line);
						line.sprintf("arcs[%d][%d].cps[%d]"
							 ".pos.y=%d\n", i, j, k,
							 (*arcs[i][j]->cps)[k]->pos.y);
						gvp.Write(line);
					}
					} else {
						line.sprintf( "arcs[%d][%d]=0\n", i, j );
						gvp.Write(line);
					}
				}
			}

			// frame separators
			line.sprintf( "numFrames=%d\n", frameEnds.size()+1 );
			gvp.Write(line);
			for (unsigned int i = 0; i < frameEnds.size(); i++) {
				line.sprintf( "frameEnds[%d].x=%d\n", i, frameEnds[i]->x );
				gvp.Write(line);
				// This can be taken from the str file.
				/*line.sprintf( "frameEnds[%d].chunkBorder=%d\n",
						  i, (frameEnds[i]->chunkBorder ? 1 : 0) );
				gvp.Write(line);*/
			}
			for (unsigned int i = 0; i < frameNameTags.size(); i++) {
				/*line.sprintf( "frames[%d].nametag.name=%s\n",
						  i, frameNameTags[i]->name.c_str() );
				gvp.Write(line);*/
				line.sprintf( "frames[%d].nametag.pos.x=%d\n",
						  i, frameNameTags[i]->pos.x );
				gvp.Write(line);
				line.sprintf( "frames[%d].nametag.pos.y=%d\n",
						  i, frameNameTags[i]->pos.y );
				gvp.Write(line);
				line.sprintf( "frames[%d].nametag.visible=%d\n",
						  i, frameNameTags[i]->isVisible() );
				gvp.Write(line);
			}

			//for all the pens, if they are not the defaults, save their colour,
			//width and style in the gvp file, style is saved as text so we have
			//to use a Map
			for(int i = 0; i < numPens; i++){
				if (!PenArray[i].pen->isDefault()){
					line.sprintf( "%s.width=%i\n", 
							PenArray[i].nameabr.c_str(),
							PenArray[i].pen->GetWidth());
					gvp.Write(line);
					line.sprintf( "%s.colour=%i,%i,%i\n", 
							PenArray[i].nameabr.c_str(),
							PenArray[i].pen->GetColour().Red(),
							PenArray[i].pen->GetColour().Green(),
							PenArray[i].pen->GetColour().Blue()
							);
					gvp.Write(line);
					//make sure we have the style in our map
					if (styleIntToStringMap.count(PenArray[i].pen->GetStyle())){
						line.sprintf( "%s.style=%s\n", 
								PenArray[i].nameabr.c_str(),
								styleIntToStringMap[PenArray[i].pen->GetStyle()].c_str());
						gvp.Write(line);
					} else {
						//XXX figure out a better way to log these errors, which shouldn't happen
						cout << "Style " << PenArray[i].pen->GetStyle() <<" is not in map and so will not be written into\n";
						cout << "the gvp File, this is a problem with the program and so you should contact the maintainer.\n";
					}
				}
			}
//			//write the font if it isn't the default font
//			if (labelFont != defaultFont) 
//				XXX: do something here

			gvpDirty = false;
			//there is no old gvp file once we have saved
			gvpFile_old = "";
		} else {
			wxString msg;
			msg.sprintf("Failed to open position file \"%s\" for writing", gvpFile.c_str());
			wxLogMessage(msg);
			//if the save wasn't sucessful then we shouldn't point to the unsucessful gvpfile
			if (gvpFile_old != "") {
				gvpFile = gvpFile_old;
				gvpFile_old = "";
			}
		}
	} else {
		SaveAs();
		//there is no old gvp file once we have saved
		gvpFile_old = "";
	}
	// We may not be dirty anymore, so update the status bar.
	onComeForward();
	//there is no old gvp file once we have saved
	gvpFile_old = "";
}

/**
 *******************************************************************
 * Prompt the user for a filename and then save the graph data to that
 * file (unless the user cancels it of course).
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The selected file may be overwritten with the current graph data.
 *
 * \note The selected file may be overwritten with the current graph data.
 *
 * \return void
 *******************************************************************/
void
StructPage::SaveAs( void )
{
	wxString defaultGvp;
	if (gvpFile == "")
		defaultGvp.sprintf("%s.gvp", strFile.c_str());
	else
		defaultGvp.sprintf("%s", gvpFile.c_str());

	wxFileDialog dlg(this, "Save position file as...",
			 defaultGvp.BeforeLast('/'), defaultGvp.AfterLast('/'),
			 "All files|*"
			 "|Position Files (*.gvp)|*.gvp"
			 "|Text Files|*.txt;*.text",
			 /* i believe that the wxCHANGE_DIR was cause problems
			  * because it would change our working directory so none
			  * of our file locations (which were all relative to the cwd)
			  * were valid
			  */
			 //wxSAVE | wxCHANGE_DIR | wxOVERWRITE_PROMPT,
			 wxSAVE | wxOVERWRITE_PROMPT,
			 wxDefaultPosition);

	dlg.SetFilterIndex(1); // show only gvps by default

	// Where do they want to save it?
	if ( dlg.ShowModal() == wxID_OK ) {
		gvpFile_old = gvpFile;
		gvpFile = dlg.GetPath();
		// save to the file
		Save();
		// put the filename wherever it needs to go
		parentFrame->SetStatusText(dlg.GetFilename(), 0);
	}
}

/**
 *******************************************************************
 * Prepare for the tab to be closed, by asking to the user to save if
 * the file is dirty. Let the caller know whether they may close the tab.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The graph data may be saved to a file.
 *
 * \note The graph data may be saved to a file.
 *
 * \return void
 *******************************************************************/
bool
StructPage::RequestClose( void )
{
	// if the file has unsaved changes
	while (gvpDirty) {
		wxString msg;
		msg.sprintf( "The file '%s' may have unsaved changes. "
				 "What would you like to do?", gvpFile.c_str() );
		wxString choices[] = { wxT("Save the changes."),
						wxT("Save the changes to a different file."),
						wxT("Discard the changes."),
						wxT("Nevermind. Don't close the file.") };
		wxSingleChoiceDialog dlg( this, msg, wxT("Save file?"), 4, choices,
					  NULL, wxOK | wxCENTRE, wxDefaultPosition );
		// prompt to save it
		dlg.ShowModal(); // no need to test for wxOK
		// if the user chooses to save it
		if (dlg.GetSelection() == 0) {
			// try saving it
			Save();
			// go back to the beginning
		} else if (dlg.GetSelection() == 1) {
			// try saving it with a different name
			SaveAs();
			// go back to the beginning
		} else if (dlg.GetSelection() == 2) {
			// else if the user discards changes, allow it to closed
			return true;
		} else if (dlg.GetSelection() == 3) {
			// else if the user cancels, do not allow it to be closed
			return false;
		}
	}
	// else allow the page to be closed
	return true;
}

/**
 *******************************************************************
 * Set all control points to have the given selected value if they are
 * endpoints attached to the node specified by the given Id.
 *
 * \param newSelected The selected value desired for the endpoints.
 * \param rvId The Id of the node whose endpoints should be (un)selected.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The control points that start an arc from or end an arc to
 * the given node will be selected (or unselected).
 *
 * \note The control points that start an arc from or end an arc to
 * the given node will be selected (or unselected).
 *
 * \return void
 *******************************************************************/
void
StructPage::setEndpointsSelected( bool newSelected, RVInfo::rvParent rvId )
{
	int n = nameVizNodeMap[rvId];
	int totalNodes = nodes.size();

	for ( int i = 0; i < totalNodes; i++ ) {
		// outgoing arcs
		if (arcs[n][i]) {
			(*arcs[n][i]->cps)[0]->setSelected(newSelected);
		}
		// incoming arcs
		if (arcs[i][n]) {
			(*arcs[i][n]->cps)[ arcs[i][n]->cps->size()-1 ]
			->setSelected(newSelected);
		}
	}
}

/**
 *******************************************************************
 * Sets the visiblity of all the arcs coming in and going out of a node
 *
 * \param newVisible The visiblity value desired for the arcs.
 * \param rvId The Id of the node whose in/out arcs should be made visible/invisible.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The arcs going in and coming out of the given node will be made visible or invisible.
 *
 * \note The arcs going in and coming out of the given node will be made visible or invisible.
 *
 * \return void
 *******************************************************************/
void
StructPage::setInOutArcsVisible( bool newVisible, RVInfo::rvParent rvId )
{
	int n = nameVizNodeMap[rvId];
	int totalNodes = nodes.size();

	for ( int i = 0; i < totalNodes; i++ ) {
		// outgoing arcs
		if (arcs[n][i]) {
			//an arc should only be made visible if both nodes are visible
			if ((newVisible && nodes[i]->isVisible()) | !newVisible){
				arcs[n][i]->setVisible(newVisible);
				for(unsigned int j = 1; j < arcs[n][i]->cps->size() - 1; j++){
					(*arcs[n][i]->cps)[j]->setVisible(newVisible);
				}
			}
		}
		// incoming arcs
		if (arcs[i][n]) {
			//an arc should only be made visible if both nodes are visible
			if ((newVisible && nodes[i]->isVisible()) | !newVisible){
				arcs[i][n]->setVisible(newVisible);
				for(unsigned int j = 1; j < arcs[i][n]->cps->size() - 1; j++){
					(*arcs[i][n]->cps)[j]->setVisible(newVisible);
				}
			}
		}
	}
}

/**
 *******************************************************************
 * Sets the highlighedness of all the arcs coming in to a node
 *
 * \param newHighlighted The highlighedness value desired for the arcs.
 * \param rvId The Id of the node whose incoming arcs should be highlighed/unhighlighed.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The arcs going in to the given node will be made highlighted or unhighlighed.
 *
 * \note The arcs going in to the given node will be made highlighted or unhighlighed.
 *
 * \return void
 *******************************************************************/
void
StructPage::setInArcsHighlighed( bool newHighlighted, RVInfo::rvParent rvId )
{
	int n = nameVizNodeMap[rvId];
	int totalNodes = nodes.size();

	for ( int i = 0; i < totalNodes; i++ ) {
		// incoming arcs
		if (arcs[i][n]) {
			(*arcs[i][n]).setHighlighted(newHighlighted);
		}
	}
}

/**
 *******************************************************************
 * Set the selected status of every item.
 *
 * \param newSelected The selected value desired for the items.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Every item will either be selected or unselected according to
 * \p newSelected.
 *
 * \note Every item will either be selected or unselected according to
 * \p newSelected.
 *
 * \return void
 *******************************************************************/
void
StructPage::setAllSelected( bool newSelected )
{
	int numNodes = nodes.size();

	// select all nodes (and associated control points)
	// (don't need to select nametags because nodes select their own nametags)
	for (int i = 0; i < numNodes; i++) {
		nodes[i]->setSelected( newSelected );
	}
	// select all other control points
	for (int i = 0; i < numNodes; i++) {
		for (int j = 0; j < numNodes; j++) {
			if (arcs[i][j]) {
				int end = arcs[i][j]->cps->size() - 1;
				assert(end > 0);
				for ( int k = 1; k < end; k++ ) {
					(*arcs[i][j]->cps)[k]->setSelected( newSelected );
				}
			}
		}
	}
	// select all frame separators
	for (unsigned int i = 0; i < frameEnds.size(); i++) {
		frameEnds[i]->setSelected(newSelected);
	}
	// select all frame nametags
	for (unsigned int i = 0; i < frameNameTags.size(); i++) {
		frameNameTags[i]->setSelected(newSelected);
	}
}

/**
 *******************************************************************
 * Toggle the selected value for all items in the given rectangle.
 *
 * \param rect The rectangle within which to toggle selected values.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The selected value for all the items in the rectangle will be
 * toggled.
 *
 * \note The selected value for all the items in the rectangle will be
 * toggled.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleSelectedInRect( const wxRect& rect )
{
	int numNodes = nodes.size();

	// select frame nametags
	for (unsigned int i = 0; i < frameNameTags.size(); i++) {
		if (frameNameTags[i]->inRect(rect))
			frameNameTags[i]->toggleSelected();
	}

	// select node nametags before nodes so that nodes will correct
	// the possibly incorrect toggle
	for (int i = 0; i < numNodes; i++) {
		if (nodeNameTags[i]->inRect(rect))
			nodeNameTags[i]->toggleSelected();
	}

	// select all nodes (and associated control points)
	for (int i = 0; i < numNodes; i++) {
		if (nodes[i]->inRect(rect))
			nodes[i]->toggleSelected();
	}
	// select all other control points
	for (int i = 0; i < numNodes; i++) {
		for (int j = 0; j < numNodes; j++) {
			if (arcs[i][j]) {
			int end = arcs[i][j]->cps->size() - 1;
			assert(end > 0);
			for ( int k = 1; k < end; k++ ) {
				if ((*arcs[i][j]->cps)[k]->inRect(rect))
				(*arcs[i][j]->cps)[k]->toggleSelected();
			}
			}
		}
	}
}

/**
 *******************************************************************
 * Don't draw any labels that are currently selected.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Currently selected labels will no longer be drawn.
 *
 * \note Currently selected labels will no longer be drawn.
 *
 * \return void
 *******************************************************************/
void
StructPage::hideSelectedLabels( void )
{
	int numNodes = nodes.size();
	for (int i = 0; i < numNodes; i++) {
		if (nodeNameTags[i]->getSelected())
			nodeNameTags[i]->setVisible(false);
	}
	for (unsigned int i = 0; i < frameNameTags.size(); i++) {
		if (frameNameTags[i]->getSelected())
			frameNameTags[i]->setVisible(false);
	}

	redraw();
	blit();
}

/**
 *******************************************************************
 * Mark all the node labels as visible and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post All of the Node Labels will be marked visible and drawn.
 *
 * \note All of the Node Labels will be marked visible and drawn.
 *
 * \return void
 *******************************************************************/
void
StructPage::showAllNodeLabels( void )
{
	for (unsigned int i = 0; i < nodes.size(); i++) {
		//only make the nametags visible if the node is visible
		if (nodes[i]->isVisible())
			nodeNameTags[i]->setVisible(true);
	}

	redraw();
	blit();
}

/**
 *******************************************************************
 * Mark all the frame labels as visible and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post All of the Frame Labels will be marked visible and drawn.
 *
 * \post All of the Frame Labels will be marked visible and drawn.
 *
 * \return void
 *******************************************************************/
void
StructPage::showAllFrameLabels( void )
{
	for (unsigned int i = 0; i < frameNameTags.size(); i++) {
		frameNameTags[i]->setVisible(true);
	}

	redraw();
	blit();
}

/**
 *******************************************************************
 * Mark all the nodes as visible and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post All of the Nodes (their labels and in/out arcs) will be marked visible and drawn.
 *
 * \note All of the Nodes (their labels and in/out arcs) will be marked visible and drawn.
 *
 * \return void
 *******************************************************************/
void
StructPage::showAllNodes( void )
{
	for (unsigned int i = 0; i < nodes.size(); i++) {
		nodes[i]->setVisible(true);
	}

	redraw();
	blit();
}

/**
 *******************************************************************
 * Toggle whether control points are drawn and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether control points will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether control points will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewCPs( void )
{
	drawCPs = !drawCPs;
	redraw();
	blit();
}

/**
 *******************************************************************
 * Toggle whether zig-zag lines are drawn and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether lines will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether lines will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewLines( void )
{
	drawLines = !drawLines;
	redraw();
	blit();
}

/**
 *******************************************************************
 * Toggle whether splines are drawn and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether splines will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether splines will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewSplines( void )
{
	drawSplines = !drawSplines;
	redraw();
	blit();
}

/**
 *******************************************************************
 * Toggle whether arrow heads are drawn and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether arrow heads will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether arrow heads will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewArrowHeads( void )
{
	drawArrowHeads = !drawArrowHeads;
	redraw();
	blit();
}

/**
 *******************************************************************
 * Toggle whether nodes are drawn and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether nodes will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether nodes will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewNodes( void )
{
	drawNodes = !drawNodes;
	redraw();
	blit();
}

/**
 *******************************************************************
 * Toggle whether grids are drawn and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether grids will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether grids will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewGrids( void )
{
	drawGrids = !drawGrids;
	redraw();
	blit();
}

/**
 *******************************************************************
 * Toggle whether direct lines are drawn and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether direct lines will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether direct lines will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewDirectLines( void )
{
	drawDirectLines = !drawDirectLines;
	redraw();
	blit();
}

/**
 *******************************************************************
 * Toggle whether frame separators are drawn and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether frame separators will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether frame separators will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewFrameSeps( void )
{
	drawFrameSeps = !drawFrameSeps;
	redraw();
	blit();
}

/**
 *******************************************************************
 * Toggle whether frame labels are drawn and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether frame labels will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether frame labels will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewFrameNames( void )
{
	drawFrameNames = !drawFrameNames;
	redraw();
	blit();
}

/**
 *******************************************************************
 * Toggle whether node labels are drawn and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether node labels will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether node labels will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewNodeNames( void )
{
	drawNodeNames = !drawNodeNames;
	redraw();
	blit();
}

/**
 *******************************************************************
 * Toggle whether tool tips are drawn and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether tool tips will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether tool tips will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \remark Tool tips in gmtkViz don't work anyway.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewToolTips( void )
{
	drawToolTips = !drawToolTips;
	wxToolTip::Enable(drawToolTips);
}
/**
 *******************************************************************
 * Toggle whether bounding box is drawn and redrawn.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether bounding box will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether bounding box will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \remark
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewBoundingBox( void )
{
	drawBoundingBox = !drawBoundingBox;
	redraw();
	blit();
}
/**
 *******************************************************************
 * Toggle whether select box is drawn and redrawn.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether select box will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \note Whether select box will be drawn or not will be toggled
 * and the screen will be refreshed.
 *
 * \remark
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleViewSelectBox( void )
{
	drawSelectBox= !drawSelectBox;
	redraw();
	blit();
}

/**
 *******************************************************************
 * Prompt the user for a new font and redraw with that font.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The font may be changed and the screen may be updated.
 *
 * \note The font may be changed and the screen may be updated.
 *
 * \return void
 *******************************************************************/
void
StructPage::changeFont( void )
{
	wxFont newFont = wxGetFontFromUser(this, labelFont);
	if (newFont.Ok()) {
		save_undo();
		labelFont = newFont;
		redraw();
		blit();
	}
}

/**
 *******************************************************************
 * Toggle whether to snap to grid when things are moved.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post Whether things are snapped to grid or not will be toggled.
 *
 * \note Whether things are snapped to grid or not will be toggled.
 *
 * \return void
 *******************************************************************/
void
StructPage::toggleSnapToGrid( void )
{
	snapToGrid = !snapToGrid;
}

/**
 *******************************************************************
 * Update the status bar with our file name and dirtiness.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The status bar will be updated.
 *
 * \note The status bar will be updated.
 *
 * \return void
 *******************************************************************/
void
StructPage::onComeForward( void )
{
	parentFrame->SetStatusText(gvpDirty?wxT("*"):wxEmptyString, 1);
	parentFrame->SetStatusText(gvpFile.length()?gvpFile:wxT("UNTITLED (and unsaved)"), 0);
}

/**
 *******************************************************************
 * Supply a reasonable name for this document.
 *
 * \param name Where to store the name.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post \p name will be updated with the appropriate string.
 *
 * \note No side effects.
 *
 * \return void
 *******************************************************************/
void
StructPage::getName( wxString& name )
{
	name = ( (gvpFile != wxEmptyString) ? gvpFile : strFile );
}

/**
 *******************************************************************
 * Update the scale at which everything is drawn/displayed.
 *
 * \param newScale The index into the gZoomMap array corresponding to
 * the desired scale.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The scale will be adjusted and the document redrawn.
 *
 * \note The scale will be adjusted and the document redrawn.
 *
 * \return void
 *******************************************************************/
void
StructPage::setScale( int newScale )
{
	displayScale = newScale;
	// The drawing buffer needs to change size.
	if (content) delete content;
	content = new wxBitmap( (int)round(canvasWidth*gZoomMap[displayScale]),
			(int)round(canvasHeight*gZoomMap[displayScale]) );
	// The scrollable area needs to change size.
	SetVirtualSize( (int)round(canvasWidth*gZoomMap[displayScale]),
			(int)round(canvasHeight*gZoomMap[displayScale]) );
	redraw();
	blit();
}

/**
 *******************************************************************
 * Prompt the user for a new canvas width and possible change it and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The canvas may be resized and redrawn.
 *
 * \note The canvas may be resized and redrawn.
 *
 * \return void
 *******************************************************************/
void
StructPage::adjustCanvasWidth( void )
{
	wxCoord width, minWidth;

	minWidth = rightMostItemX() + NODE_RADIUS;

	width = wxGetNumberFromUser( wxT("How wide should the canvas be?"),
			wxT("Width: "), wxT("Change Canvas Width"),
			getWidth(), minWidth, getWidth()*2 );

	if (width != -1) {
		canvasWidth = width;
		setScale(getScale());
	}
}

/**
 *******************************************************************
 * Prompt the user for a new canvas height and possible change it and redraw.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The canvas may be resized and redrawn.
 *
 * \note The canvas may be resized and redrawn.
 *
 * \return void
 *******************************************************************/
void
StructPage::adjustCanvasHeight( void )
{
	wxCoord height, minHeight;

	minHeight = bottomMostItemY() + NODE_RADIUS;

	height = wxGetNumberFromUser( wxT("How tall should the canvas be?"),
				  wxT("Height: "), wxT("Change Canvas Height"),
				 getHeight(), minHeight, getHeight()*2 );

	if (height != -1) {
		canvasHeight = height;
		setScale(getScale());
	}
}

/**
 *******************************************************************
 * Return the X coordinate of the right-most item.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post No real changes.
 *
 * \note No real changes.
 *
 * \return the X coordinate of the right-most item
 *******************************************************************/
long
StructPage::rightMostItemX( void )
{
	long xMax = 0;
	int numNodes = nodes.size();
	// check all nodes (and associated control points)
	for (int i = 0; i < numNodes; i++) {
		xMax = ( nodes[i]->center.x <= xMax ? xMax : nodes[i]->center.x );
	}
	// check all node nametags
	for (int i = 0; i < numNodes; i++) {
		xMax = ( nodeNameTags[i]->pos.x <= xMax
			 ? xMax
			 : nodeNameTags[i]->pos.x );
	}
	// check all other control points
	for (int i = 0; i < numNodes; i++) {
		for (int j = 0; j < numNodes; j++) {
			if (arcs[i][j]) {
				int end = arcs[i][j]->cps->size() - 1;
				assert(end > 0);
				for ( int k = 1; k < end; k++ ) {
					( (*arcs[i][j]->cps)[k]->pos.x <= xMax
					  ? xMax
					  : (*arcs[i][j]->cps)[k]->pos.x );
				}
			}
		}
	}
	// check all frame separators
	for (unsigned int i = 0; i < frameEnds.size(); i++) {
		xMax = ( frameEnds[i]->x <= xMax ? xMax : frameEnds[i]->x );
	}
	// check all frame nametags
	for (unsigned int i = 0; i < frameNameTags.size(); i++) {
		xMax = ( frameNameTags[i]->pos.x <= xMax
			 ? xMax
			 : frameNameTags[i]->pos.x );
	}

	return xMax;
}

/**
 *******************************************************************
 * Return the Y coordinate of the bottom-most item.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post No real changes.
 *
 * \note No real changes.
 *
 * \return the Y coordinate of the bottom-most item
 *******************************************************************/
long
StructPage::bottomMostItemY( void )
{
	long yMax = 0;
	int numNodes = nodes.size();
	// check all nodes (and associated control points)
	for (int i = 0; i < numNodes; i++) {
		yMax = ( nodes[i]->center.y <= yMax ? yMax : nodes[i]->center.y );
	}
	// check all node nametags
	for (int i = 0; i < numNodes; i++) {
		yMax = ( nodeNameTags[i]->pos.y <= yMax
			 ? yMax
			 : nodeNameTags[i]->pos.y );
	}
	// check all other control points
	for (int i = 0; i < numNodes; i++) {
		for (int j = 0; j < numNodes; j++) {
			if (arcs[i][j]) {
				int end = arcs[i][j]->cps->size() - 1;
				assert(end > 0);
				for ( int k = 1; k < end; k++ ) {
					( (*arcs[i][j]->cps)[k]->pos.y <= yMax
					  ? yMax
					  : (*arcs[i][j]->cps)[k]->pos.y );
				}
			}
		}
	}
	// check all frame nametags
	for (unsigned int i = 0; i < frameNameTags.size(); i++) {
		yMax = ( frameNameTags[i]->pos.y <= yMax
			 ? yMax
			 : frameNameTags[i]->pos.y );
	}

	return yMax;
}

/**
 *******************************************************************
 * Save the state of the StructPage so that it can be restored later.
 *
 * \pre The StructPage should be fully initialized.
 *
 * \post The undo_stack will have a new item on it containing the current
 * state.
 *
 * \note The undo_stack will have a new item on it containing the current
 * state.
 *
 * \return nothing
 *******************************************************************/
void
StructPage::save_state(std::stack<StructPage_state *> &state_stack)
{
	//make a new state structure
	StructPage_state * state = new StructPage_state;
	//save the state
	for (int i = 0; i < numPens; i++){
		state->PenArray[i].pen = new gmtkPen(*(PenArray[i].pen));
		state->PenArray[i].nameabr = PenArray[i].nameabr;
		state->PenArray[i].namelong = PenArray[i].namelong;
	}
	state->labelFont = labelFont;
	state->drawSelectBox = drawSelectBox;
	state->drawCPs = drawCPs;
	state->drawLines = drawLines;
	state->drawSplines = drawSplines;
	state->drawArrowHeads = drawArrowHeads;
	state->drawNodes = drawNodes;
	state->drawGrids = drawGrids;
	state->drawDirectLines = drawDirectLines;
	state->drawFrameSeps = drawFrameSeps;
	state->drawNodeNames = drawNodeNames;
	state->drawFrameNames = drawFrameNames;
	state->drawToolTips = drawToolTips;
	state->drawBoundingBox = drawBoundingBox;
	state->displayScale = displayScale;
	state->canvasWidth = canvasWidth;
	state->canvasHeight = canvasHeight;
	state->snapToGrid = snapToGrid;

	std::vector< VizArc* > curVec;

	int numNodes = nodes.size();
	state->nodes.clear();
	state->arcs.clear();
	for (int i = 0; i < numNodes; i++){
		curVec.clear();
		state->arcs.push_back(curVec);
		state->nodes.push_back(new VizNode(*nodes[i]));
		state->nodeNameTags.push_back(&(state->nodes[i]->nametag));
		for (int j = 0; j < numNodes; j++){
			if (arcs[i][j] != NULL)
				state->arcs[i].push_back(new VizArc(*(arcs[i][j])));
			else
				state->arcs[i].push_back(NULL);
		}
	}

	//just incase
	state->frameEnds.clear();
	state->frameNameTags.clear();
	for (unsigned int i = 0; i < frameEnds.size(); i++)
		state->frameEnds.push_back(new VizSep(*frameEnds[i]));
	for (unsigned int i = 0; i < frameNameTags.size(); i++)
		state->frameNameTags.push_back(new NameTag(*frameNameTags[i]));

	//push the state onto the given stack
	state_stack.push(state);
}

bool
StructPage::restore_state(std::stack<StructPage_state *> &state_stack)
{
	if(state_stack.empty()){
		return false;	//we cannot restore state
	}
	StructPage_state * state;
	state = state_stack.top();
	for (int i = 0; i < numPens; i++){
		PenArray[i].pen->SetColour(state->PenArray[i].pen->GetColour());
		PenArray[i].pen->SetWidth(state->PenArray[i].pen->GetWidth());
		PenArray[i].pen->SetStyle(state->PenArray[i].pen->GetStyle());
		delete state->PenArray[i].pen;
		state->PenArray[i].pen = NULL;
		//XXX do we need to do this below?
//		delete state->PenArray[i].nameabr;
//		delete state->PenArray[i].namelong;
	}
	labelFont = state->labelFont;
	drawSelectBox = state->drawSelectBox;
	drawCPs = state->drawCPs;
	drawLines = state->drawLines;
	drawSplines = state->drawSplines;
	drawArrowHeads = state->drawArrowHeads;
	drawNodes = state->drawNodes;
	drawGrids = state->drawGrids;
	drawDirectLines = state->drawDirectLines;
	drawFrameSeps = state->drawFrameSeps;
	drawNodeNames = state->drawNodeNames;
	drawFrameNames = state->drawFrameNames;
	drawToolTips = state->drawToolTips;
	drawBoundingBox = state->drawBoundingBox;
	displayScale = state->displayScale;
	canvasWidth = state->canvasWidth;
	canvasHeight = state->canvasHeight;
	snapToGrid = state->snapToGrid;

	nodes = state->nodes;
	nodeNameTags = state->nodeNameTags;
	arcs = state->arcs;
	frameNameTags = state->frameNameTags;
	frameEnds = state->frameEnds;
	state_stack.pop();
	delete state;
	return true;
}

void
StructPage::save_undo()
{ 
	save_state(undo_stack);
	//if we save an undo state then we have no need for the
	//redo stack
	while(!redo_stack.empty()){
		delete redo_stack.top();
		redo_stack.pop();
	}
}

bool
StructPage::undo()
{ 
	if(undo_stack.empty())
		return false;
	save_state(redo_stack);
	return restore_state(undo_stack);
}

bool
StructPage::redo()
{ 
	if(redo_stack.empty())
		return false;
	save_state(undo_stack);
	return restore_state(redo_stack);
}

/**
 *******************************************************************
 * Construct a NameTag.
 *
 * \param newPos The location of the upper right corner of the label.
 * \param newName The text to be displayed in the label.
 *
 * \pre Valid parameters.
 *
 * \post The NameTag should exist.
 *
 * \note No side effects.
 *
 * \return Nothing.
 *******************************************************************/
NameTag::NameTag( const wxPoint& newPos, const wxString& newName )
	: pos(newPos), name(newName)
{
	visible = true;
	selected = false;
	frame = -1; //indicates that we haven't given a frame value
}
/**
 *******************************************************************
 * Construct a NameTag from an existing NameTag.
 *
 * \param original, the nametag to copy
 *
 * \pre Valid NameTag.
 *
 * \post The NameTag should exist, it will be a copy of original.
 *
 * \note No side effects.
 *
 * \return a NameTag object.
 *******************************************************************/
NameTag::NameTag( const NameTag &original)
{
	visible = original.visible;
	pos = original.pos;
	name = original.name;
	selected = false;
	frame = original.frame;
}

/**
 *******************************************************************
 * Drawing routine for the NameTag.
 *
 * \param dc A pointer to the device context in which to draw the NameTag.
 *
 * \pre \p dc should be ready.
 *
 * \post \p dc will be updated.
 *
 * \note \p dc will be updated.
 *
 * \return void
 *******************************************************************/
void
NameTag::draw( wxDC *dc )
{
	if (!visible) 
		return;
	dc->GetTextExtent(name, &size.x, &size.y);
	if (getSelected())
		dc->DrawRectangle(pos.x, pos.y, size.x, size.y);
	dc->DrawText(name, pos);
}

/**
 *******************************************************************
 * Determine whether a given point lies on the NameTag.
 *
 * \param pt The point that might be on the NameTag.
 *
 * \pre The NameTag should have been drawn at least once (no we know
 * how big it is in the device context).
 *
 * \post No real changes.
 *
 * \note No side effects.
 *
 * \return \c true if the point lies on the the NameTag, \c false if
 * it doesn't.
 *******************************************************************/
bool
NameTag::onMe( const wxPoint& pt )
{
	wxRect temp(pos.x, pos.y, size.x, size.y);
	return temp.Inside(pt);
}

/**
 *******************************************************************
 * Determine whether the NameTag lies within the given rectangle.
 *
 * \param rect The rectangle which might contain the NameTag.
 *
 * \pre \p rect and the NameTag should be valid.
 *
 * \post No real changes.
 *
 * \note No side effects.
 *
 * \return \c true if the center of the NameTag is inside the
 * rectangle, or false otherwise.
 *******************************************************************/
bool
NameTag::inRect( const wxRect& rect )
{
	return rect.Inside(pos.x + size.x/2, pos.y + size.y/2);
}

/**
 *******************************************************************
 * Constructs a VizNode.
 *
 * \param pos The center of the node.
 * \param newRvi The RVInfo (name and frame) for this node.
 * \param newPage The StructPage owning this node.
 *
 * \pre Valid parameters.
 *
 * \post This VizNode will be initialized along with its NameTag.
 *
 * \note No side effects.
 *
 * \return Nothing.
 *******************************************************************/

VizNode::VizNode( const wxPoint& pos, RVInfo *newRvi, StructPage *newPage )
	: nametag(pos, wxT(newRvi->name.c_str()))
{
	center.x = pos.x;
	center.y = pos.y;
	//copy the newRvi object (deep copy)
	rvi = new RVInfo(*newRvi);
	rvId.first = rvi->name;
	rvId.second = rvi->frame;
	page = newPage;
	/* This may provide tooltips on Windows, but not GTK because the
	 * window won't actually be transparent, and takes events. */
	tipWin = new wxWindow( newPage, -1, wxPoint(pos.x - NODE_RADIUS,
						pos.y - NODE_RADIUS),
			   wxSize(2*NODE_RADIUS, 2*NODE_RADIUS),
			   wxTRANSPARENT_WINDOW );
	tipWin->Hide();
	tipWin->SetToolTip(wxT(rvId.first.c_str()));
	visible = true;
	highlight_state = off;
	nametag.SetFrame(rvi->frame);
}

/**
 *******************************************************************
 * Constructs a VizNode, which is a copy of a given VizNode
 *
 * \param original, a VizNode to be copied.
 *
 * \pre Valid parameters.
 *
 * \post This VizNode will be created, it will have a Valid NameTag
 * and they will be copies of the original VizNode
 *
 * \note No side effects.
 *
 * \return Nothing.
 *******************************************************************/

VizNode::VizNode(const VizNode &original) : nametag(original.nametag)
{
	center.x = original.center.x;
	center.y = original.center.y;
	rvi = new RVInfo(*(original.rvi));
	rvId.first = original.rvId.first;
	rvId.second = original.rvId.second;
	page = original.page;
	tipWin = original.tipWin;
	visible = original.visible;
	highlight_state = off;
}


/**
 *******************************************************************
 * Destroy the VizNode in a spectacular no-op.
 *
 * \pre The VizNode should be fully initialized.
 *
 * \post The VizNode will no longer be valid. It's memory will be freed.
 *
 * \note The VizNode will no longer be valid. It's memory will be freed.
 *
 * \return Nothing.
 *******************************************************************/
VizNode::~VizNode( void )
{
	if (rvi != NULL)
		delete rvi;
	rvi = NULL;
	// don't need to destroy tipWin because it's parent will
}

/**
 *******************************************************************
 * Draw the VizNode to the given device context.
 *
 * \param dc The device context in which to draw.
 *
 * \pre The device context must be valid.
 *
 * \post The device context will be updated.
 *
 * \note The device context will be updated.
 *
 * \return void
 *******************************************************************/
void
VizNode::draw( wxDC *dc )
{
	if (!visible)
		return;
	wxBrush oldBrush = dc->GetBrush();
	if (rvi->rvDisp == RVInfo::d_observed) {
		dc->SetBrush(*wxLIGHT_GREY_BRUSH);
	}

	if (highlight_state != off){
		wxPen oldPen = dc->GetPen();
		dc->SetPen(page->highlightPen);
		dc->DrawCircle(center, NODE_RADIUS);
		dc->SetPen(oldPen);
	} else
		dc->DrawCircle(center, NODE_RADIUS);

	dc->SetBrush(oldBrush);
	if (getSelected()) {
		dc->DrawRectangle( center.x-NODE_RADIUS, center.y-NODE_RADIUS,
				3*ACTUAL_SCALE, 3*ACTUAL_SCALE );
		dc->DrawRectangle( center.x+NODE_RADIUS, center.y-NODE_RADIUS,
				-3*ACTUAL_SCALE, 3*ACTUAL_SCALE );
		dc->DrawRectangle( center.x-NODE_RADIUS, center.y+NODE_RADIUS,
				3*ACTUAL_SCALE, -3*ACTUAL_SCALE );
		dc->DrawRectangle( center.x+NODE_RADIUS, center.y+NODE_RADIUS,
				-3*ACTUAL_SCALE, -3*ACTUAL_SCALE );
	}
}

/**
 *******************************************************************
 * Determine whether a point lies on this VizNode.
 *
 * \param pt The point which might be on this VizNode.
 *
 * \pre The point must be valid.
 *
 * \post No real changes.
 *
 * \note No side effects.
 *
 * \return \c true if the point is on the VizNode, or \c false otherwise.
 *******************************************************************/
bool
VizNode::onMe( const wxPoint& pt )
{
	int dx = (pt.x - center.x), dy = (pt.y - center.y);
	return ( hypot(dx, dy) <= NODE_RADIUS );
}

/**
 *******************************************************************
 * Set the selected state of this node and any arc endpoints leading
 * to or from it.
 *
 * \param newSelected The new selected value for this node (and
 * endpoints).
 *
 * \pre The StructPage and VizNode should be fully initialized.
 *
 * \post The node and enpoints will have the specified selected value.
 *
 * \note The node and enpoints will have the specified selected value.
 *
 * \return void
 *******************************************************************/
void
VizNode::setSelected( bool newSelected ) {
	//if it isn't visible we shouldn't be able to select it
	if (!visible){
		newSelected = false;
	}
	selected = newSelected;
	// the selected state of the node must be the same as all the
	// control points belonging to the node
	page->setEndpointsSelected( newSelected, rvId );
	// select the node's nametag as well
	nametag.setSelected( newSelected );
}

/**
 *******************************************************************
 * Set the visibility state of this node, it's nametag and any arc endpoints
 * leading to or from it.
 *
 * \param newVisible The new visiblity value for this node, nametag (and
 * endpoints).
 *
 * \pre The StructPage and VizNode should be fully initialized.
 *
 * \post The node and enpoints will have the specified selected value.
 *
 * \note The node and enpoints will have the specified selected value.
 *
 * \return void
 *******************************************************************/
void
VizNode::setVisible( bool newVisible ) {
	visible = newVisible;
	//the visiblity of the arcs going in and out should be the same as
	//the node
	page->setInOutArcsVisible( newVisible, rvId );
	// change the visiblity of the node's nametag as well
	nametag.setVisible( newVisible );
}

/**
 *******************************************************************
 * Steps to the next highlight state off-> on, on -> off, children -> parents
 * parents -> both, both->off
 *
 * \param void
 *
 * \pre VizNode should be fully initialized.
 *
 * \post No side effects.
 *
 * \note No side effects.
 *
 * \return the next highlight state
 *******************************************************************/
highlight_states 
VizNode::getNextHighlightState(){
	switch (highlight_state) {
		case off: return highlight_state = on;
		case on: return highlight_state = off;
		case children: return highlight_state = parents;
		case parents: return highlight_state = both;
		case both: return highlight_state = off;
		default: return highlight_state = off;
	}
}


/**
 *******************************************************************
 * Determine whether a point lies on this ControlPoint
 *
 * \param pt The point that might be on this ControlPoint.
 *
 * \pre The point must be valid.
 *
 * \post No real changes.
 *
 * \note No side effects.
 *
 * \return \c true if the point is on the ControlPoint, or false otherwise.
 *******************************************************************/
bool
ControlPoint::onMe( const wxPoint& pt )
{
	return pos.x-1*ACTUAL_SCALE <= pt.x && pt.x <= pos.x+1*ACTUAL_SCALE
		&& pos.y-1*ACTUAL_SCALE <= pt.y && pt.y <= pos.y+1*ACTUAL_SCALE;
}

/**
 *******************************************************************
 * Construct the ControlPoint.
 *
 * \param pt The point on which to put this ControlPoint.
 *
 * \pre The point should be valid.
 *
 * \post The ControlPoint should be valid except that its arc pointer
 * does not yet point to its parent. This must be set manually.
 *
 * \note No side effects.
 *
 * \return Nothing.
 *******************************************************************/
ControlPoint::ControlPoint( const wxPoint& pt ) : Selectable()
{
	pos.x = pt.x;
	pos.y = pt.y;
	arc = NULL;
}

/**
 *******************************************************************
 * Construct a VizArc.
 *
 * \param newCps A vector of ControlPoints for the arc's spline to
 * follow.
 * \param newPage The StructPage owning this VizArc.
 *
 * \pre Valid parameters. The StructPage should be fully initialized.
 *
 * \post The VizArc should be fully initialized. It will adopt its
 * ControlPoints, pointing their arc pointers to itself. It will also
 * create the wxList of wxPoints. If you add or remove ControlPoints
 * later, the wxList and arc pointers will have to be updated manually.
 *
 * \note No side effects.
 *
 * \return Nothing.
 *******************************************************************/
VizArc::VizArc( std::vector< ControlPoint* > *newCps, StructPage *newPage )
{
	cps = newCps;
	page = newPage;

	comb_type = undefined;

	for (unsigned int i = 0; i < cps->size(); i++) {
		(*cps)[i]->arc = this;
	}
	
	int numPoints = cps->size();
	points = new wxList;
	for (int i = 0; i < numPoints; i++) {
		points->Append( (wxObject*)&(*cps)[i]->pos );
	}
	visible = true;
	highlighted = false;
}

VizArc::VizArc(const VizArc &original)
{
	page = original.page;

	comb_type = original.comb_type;

	visible = original.visible;
	highlighted = false;

	cps = new std::vector< ControlPoint* >;

	(*cps).clear();
	for (unsigned int i = 0; i < original.cps->size(); i++) {
		cps->push_back(new ControlPoint(*(*(original.cps))[i]));
		(*cps)[i]->arc = this;
		(*cps)[i]->pos.x = (*(original.cps))[i]->pos.x;
		(*cps)[i]->pos.y = (*(original.cps))[i]->pos.y;
	}
	int numPoints = cps->size();
	points = new wxList;
	for (int i = 0; i < numPoints; i++) {
		points->Append( (wxObject*)&(*cps)[i]->pos );
	}
}

/**
 *******************************************************************
 * Destroy the VizArc, freeing dynamically allocated data.
 *
 * \pre The VizArc should be fully initialized.
 *
 * \post The ControlPoints will be deleted as well as the wxList and vector.
 *
 * \note The ControlPoints will be deleted as well as the wxList and vector.
 *
 * \return Nothing.
 *******************************************************************/
VizArc::~VizArc( void )
{
	for (int k = cps->size() - 1; k >= 0; k--) {
		delete (*cps)[k];
		(*cps)[k] = NULL;
	}
	delete points;
	points = NULL;
	delete cps;
	cps = NULL;
}

/**
 *******************************************************************
 * Draw this arc to the given device context.
 *
 * \param dc The device context on which to draw.
 *
 * \pre The StructPage and wxDC should be fully initialized.
 *
 * \post The device context will be updated.
 *
 * \note The device context will be updated.
 *
 * \return void
 *******************************************************************/
void
VizArc::draw( wxDC *dc, int drawFlags )
{
	if (!visible)
		return;

	wxPen oldPen = dc->GetPen();
	if (drawFlags & DRAW_ARCS) {
		wxPen *pen = NULL;

		/* Set the pen appropriately */
		if (highlighted)
			pen = &page->highlightPen;
		else {
			switch (comb_type){
				case determin:
					pen = &page->detPen;
					break;
				case random:
					pen = &page->randPen;
					break;
				case det_rand:
					pen = &page->det_randPen;
					break;
				case switching:
					pen = &page->switchingPen;
					break;
				case det_switch:
					pen = &page->det_switchPen;
					break;
				case rand_switch:
					pen = &page->rand_switchPen;
					break;
				case det_rand_switch:
					pen = &page->det_rand_switchPen;
					break;
				default:
					cout << "error, pen has no type\n";
					return;
			}
		}

		dc->SetPen(*pen);

		/* Draw whatever splines, lines, or direct lines the StructPage
		 * calls for. */
		if ( page->getViewSplines() )
			dc->DrawSpline(points);
		if ( page->getViewLines() )
			dc->DrawLines(points);
		if ( page->getViewDirectLines() )
			dc->DrawLine((*cps)[0]->pos, (*cps)[cps->size()-1]->pos);

		if ( page->getViewArrowHeads() ) {
			// draw the arrow (quite a complicated process really)
			wxPoint arrow[3];
			int opp, adj;
			double hyp;
			opp = (*cps)[cps->size()-1]->pos.y - (*cps)[cps->size()-2]->pos.y;
			adj = (*cps)[cps->size()-1]->pos.x - (*cps)[cps->size()-2]->pos.x;
			hyp = hypot( opp, adj );

			arrow[0].x = (*cps)[cps->size()-1]->pos.x -
				(int)(NODE_RADIUS * adj / hyp);
			arrow[0].y = (*cps)[cps->size()-1]->pos.y -
				(int)(NODE_RADIUS * opp / hyp);

			arrow[1].x = arrow[0].x - (int)round(ARROW_LEN * adj / hyp) +
				(int)round(ARROW_WID * opp / hyp);
			arrow[1].y = arrow[0].y - (int)round(ARROW_LEN * opp / hyp) -
				(int)round(ARROW_WID * adj / hyp);
			arrow[2].x = arrow[0].x - (int)round(ARROW_LEN * adj / hyp) -
				(int)round(ARROW_WID * opp / hyp);
			arrow[2].y = arrow[0].y - (int)round(ARROW_LEN * opp / hyp) +
				(int)round(ARROW_WID * adj / hyp);

			wxBrush oldBrush = dc->GetBrush();
			wxBrush newBrush(oldBrush);
			newBrush.SetColour(pen->GetColour());
			dc->SetBrush(newBrush);
			dc->DrawPolygon( 3, arrow );
			dc->SetBrush(oldBrush);
		}
		dc->SetPen(oldPen);
	}

	if (drawFlags & DRAW_CPS) {
		// draw the control points
		if (highlighted)
			dc->SetPen(page->highlightPen);
		else
			dc->SetPen(page->controlPointPen);

		int end = points->GetCount() - 1;
		assert(end > 0);
		wxNode *node = points->GetFirst();
		node = node->GetNext();
		for ( int i = 1; i < end; i++, node = node->GetNext() ) {
			assert(node != NULL);
			wxPoint *p = (wxPoint *)node->GetData();
			if ( (*cps)[i]->getSelected() )
				dc->DrawRectangle( p->x-(3*ACTUAL_SCALE)/2,
						p->y-(3*ACTUAL_SCALE)/2,
						3*ACTUAL_SCALE, 3*ACTUAL_SCALE );
			if ( page->getViewCPs() ) {
				dc->DrawRectangle( p->x-ACTUAL_SCALE, p->y-ACTUAL_SCALE,
						2*ACTUAL_SCALE, 2*ACTUAL_SCALE );
			}
		}
		dc->SetPen(oldPen);
	}
}


/**
 *******************************************************************
 * Construct the VizSep.
 *
 * \param newX The x coordinate of the desired frame separator.
 * \param newPage The StructPage owning this VizSep.
 * \param newChunkBorder Is this frame separator at either end of the
 * chunk?
 *
 * \pre Valid parameters. The StructPage should be fully initialized.
 *
 * \post The VizSep should be valid.
 *
 * \note No side effects.
 *
 * \return Nothing.
 *******************************************************************/
VizSep::VizSep( wxCoord newX, StructPage *newPage, FrameSepType newSepType )
{
	x = newX;
	page = newPage;
	sepType = newSepType;
}
/**
 *******************************************************************
 * Construct the VizSep, a copy of a give VizSep
 *
 * \param a VizSep to copy
 *
 * \pre Valid parameters.
 *
 * \post The VizSep should be a valid copy of original.
 *
 * \note No side effects.
 *
 * \return Nothing.
 *******************************************************************/

VizSep::VizSep(const VizSep &original)
{
	x = original.x;
	page = original.page;
	sepType = original.sepType;
}

/**
 *******************************************************************
 * Draw the VizSep on the given device context.
 *
 * \param dc The device context on which to draw.
 *
 * \pre The device context should be fully initialized.
 *
 * \post The device context will be updated.
 *
 * \note The device context will be updated.
 *
 * \return void
 *******************************************************************/
void
VizSep::draw( wxDC * dc )
{
	// pick a pen
	wxPen *pen = &page->frameBorderPen;
	if (sepType == BEGIN_CHUNK || sepType == END_CHUNK)
		pen = &page->chunkBorderPen;

	// remember the old pen
	wxPen oldPen = dc->GetPen();
	dc->SetPen(*pen);
	// draw the frame separator
	dc->DrawLine( x, 0, x, page->getHeight() );
	dc->SetPen(oldPen);

	// maybe draw handles
	if (getSelected()) {
		dc->DrawRectangle( x-1*ACTUAL_SCALE, -1*ACTUAL_SCALE,
					3*ACTUAL_SCALE, 3*ACTUAL_SCALE );
		dc->DrawRectangle( x-1*ACTUAL_SCALE, page->getHeight()-2*ACTUAL_SCALE,
					3*ACTUAL_SCALE, 3*ACTUAL_SCALE );
	}
}

/**
 *******************************************************************
 * Determine whether the given point is on this VizSep.
 *
 * \param pt The point which might be on this VizSep.
 *
 * \pre The point should be valid.
 *
 * \post No real changes.
 *
 * \note No side effects.
 *
 * \return \c true if the point is on the VizSep, or \c false otherwise.
 *******************************************************************/
bool
VizSep::onMe( const wxPoint& pt )
{
	return x-1*ACTUAL_SCALE <= pt.x && pt.x <= x+1*ACTUAL_SCALE
	&& 0 <= pt.y && pt.y <= page->getHeight();
}


/**
 *******************************************************************
 * Construct the GmtkPrintout.
 *
 * \param newPage The StructPage whose graph should be drawn.
 * \param title The name of this printout.
 *
 * \pre The StructPage and title should be fully initialized.
 *
 * \post The GmtkPrintout should be valid.
 *
 * \note Use of this object may cause the StructPage to be drawn with
 * a different scale to a different device context, which in turn may
 * alter the size of NameTag objects until it is drawn to its usual
 * device context with its usual scale again.
 *
 * \return Nothing.
 *******************************************************************/
GmtkPrintout::GmtkPrintout(StructPage * newPage, const wxChar *title)
	: wxPrintout(title)
{
	page = newPage;
}

/**
 *******************************************************************
 * Print the requested page to whatever device context is given by
 * GetDC().
 *
 * \param page Should be 1, because we only have one page.
 *
 * \pre The StructPage and GmtkPrintout should be fully initialized.
 *
 * \post The device context will be updated and the NameTag objects
 * may have changed size according to this device context and
 * scale. To revert their sizes, draw to the usual device context with
 * the usual scale.
 *
 * \note The device context will be updated and the NameTag objects
 * may have changed size according to this device context and
 * scale. To revert their sizes, draw to the usual device context with
 * the usual scale.
 *
 * \return \c true on success, or \c false on failure.
 *******************************************************************/
bool
GmtkPrintout::OnPrintPage(int page)
{
	wxDC *dc = GetDC();
	if (dc) {
		if (page == 1)
			DrawPageOne(dc);
		else 
			return false;
		
		dc->SetDeviceOrigin(0, 0);
		dc->SetUserScale(1.0, 1.0);
		
		//wxChar buf[200];
		//wxSprintf(buf, wxT("PAGE %d"), page);
		// dc->DrawText(buf, 10, 10);
		
		return true;
	}
	else
		return false;
}

/**
 *******************************************************************
 * Disclose the page ranges of this document (which only has one page).
 *
 * \pre The GmtkPrintout should be fully initialized.
 *
 * \post The parameters will be filled with 1s.
 *
 * \note No side effects.
 *
 * \return void
 *******************************************************************/
void
GmtkPrintout::GetPageInfo( int *minPage, int *maxPage,
			   int *selPageFrom, int *selPageTo )
{
	*minPage = 1;
	*maxPage = 1;
	*selPageFrom = 1;
	*selPageTo = 1;
}

/**
 *******************************************************************
 * Draw the graph to the given device context, scaling it to fit in
 * the drawable area.
 *
 * \param dc The device context on which to draw.
 *
 * \pre The GmtkPrintout and wxDC should be fully initialized.
 *
 * \post The device context will be updated and the NameTag objects
 * may have changed size according to this device context and
 * scale. To revert their sizes, draw to the usual device context with
 * the usual scale.
 *
 * \note The device context will be updated and the NameTag objects
 * may have changed size according to this device context and
 * scale. To revert their sizes, draw to the usual device context with
 * the usual scale.
 *
 * \return void
 *******************************************************************/
void
GmtkPrintout::DrawPageOne(wxDC *dc)
{
/* You might use THIS code if you were scaling
* graphics of known size to fit on the page.
	*/
	int w, h;
	
	float maxX = page->getWidth();
	float maxY = page->getHeight();
	
	// Let's have at least 10 device units margin
	float marginX = 10;
	float marginY = 10;
	
	// Add the margin to the graphic size
	maxX += (2*marginX);
	maxY += (2*marginY);
	
	// Get the size of the DC in pixels
	dc->GetSize(&w, &h);
	
	// Calculate a suitable scaling factor
	float scaleX=(float)(w/maxX);
	float scaleY=(float)(h/maxY);
	
	// Use x or y scaling factor, whichever fits on the DC
	float actualScale = wxMin(scaleX,scaleY);
	
	// Calculate the position on the DC for centring the graphic
	float posX = (float)((w - (page->getWidth()*actualScale))/2.0);
	float posY = (float)((h - (page->getHeight()*actualScale))/2.0);
	
	// Set the scale and origin
	dc->SetUserScale(actualScale, actualScale);
	dc->SetDeviceOrigin( (long)posX, (long)posY );
	//dc->SetUserScale(1.0, 1.0);

	//wxFont oldFont = dc->GetFont();
	//wxFont scaledFont(oldFont);
	//scaledFont.SetPointSize((int)(actualScale*oldFont.GetPointSize()));
	//dc->SetFont(scaledFont);

	//turn off controlpts and selectBox if they're viewable
	bool had_view_cps = page->getViewCPs();
	if(had_view_cps)
		page->toggleViewCPs();
	bool had_view_select_box = page->getViewSelectBox();
	if(had_view_select_box)
		page->toggleViewSelectBox();

	//deselect everything
	page->setAllSelected(false);

	page->draw(*dc);

	//dc->SetFont(oldFont);
	
	//revert to the previous settings with the controlpts and selectBox
	if(had_view_cps)
		page->toggleViewCPs();
	if(had_view_select_box)
		page->toggleViewSelectBox();
}

 /*******************************************************************
  * Create the Contents of the Help window
 *
 * \param 
 *
 * \pre The GFrame should be fully initialized
 *
 * \post The GmktHelp object will contain data and a button
 *
 * \note The GmktHelp object will contain data and a button
 * 
 * \return void
 *******************************************************************/

void
GmtkHelp::doLayout()
{
	wxBoxSizer* help_sizer = new wxBoxSizer(wxVERTICAL);
	wxTextCtrl* help_msg = new wxTextCtrl(this, -1, "", wxDefaultPosition,
			wxSize(800,300), wxTE_READONLY | wxTE_MULTILINE | wxTE_LEFT | wxTE_DONTWRAP);

	wxTextAttr normal(*wxBLACK,wxNullColour,wxFont(12, wxROMAN, wxNORMAL, wxNORMAL));
	wxTextAttr title(*wxBLACK,wxNullColour,wxFont(14, wxROMAN, wxNORMAL, wxBOLD));
	wxTextAttr bold(*wxBLACK,wxNullColour,wxFont(12, wxROMAN, wxNORMAL, wxBOLD));
	wxTextAttr italic(*wxBLACK,wxNullColour,wxFont(12, wxROMAN, wxITALIC, wxNORMAL));

	help_msg->SetDefaultStyle(title);
	help_msg->AppendText("Gmtkviz Help Sheet:\n\n");

	help_msg->SetDefaultStyle(bold);
	help_msg->AppendText("Keyboard Commands:\n");

	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'a'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : If the mouse pointer is over a node 'a' will pop up node info\n");

	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t's'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : If the mouse pointer is over a node or frame label 's' will toggle the node or "
		"frame label's visibility\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'w'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : toggles the visibility of node and frame labels in the current selection\n");
	
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'd'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : If the mouse pointer is over a node will toggle that node, its label and its edge's visibility\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'e'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : If there is an active selection will toggle the visibility of the selected nodes and their labels and edges\n");

	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'f'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : will step through highlighting modes if the mouse pointer is over a visible node.\n"
			"\t\t\tFirst: the node, its outgoing edges and its children will be highlighted.\n"
			"\t\t\tSecond: the node, its incoming edges and its parents will be highlighted.\n"
			"\t\t\tThird: the node, its edges, its children and its parents will be highlighted.\n"
			"\t\t\tFourth: all highlighting will be turned off\n"
			"\t\tturns off all highlighting if the mouse pointer is not on a visible node\n");

	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'z'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : selects all the nodes in the frame under the mouse pointer.\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'Z'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : adds the nodes in the frame under the mouse pointer to the current selection\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'x'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : invert selection.  Everything that was selected is deselected and everything that was not selected is selected.\n");

	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'delete'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : deletes selected contropoint(s). NOTE: this will not delete anything except control points\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'r'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : deletes selected contropoint(s). NOTE: this will not delete anything except control points\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'g'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : undo... steps through the undo stack\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'G'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : redo... steps through the redo stack, only exits after an undo and before any other undoable actions\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'ctrl-a'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : selects everything that is selectable.\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'ctrl-d'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : makes all nodes visible.\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'ctrl-e'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : brings up the \"Print to EPS file\" dialog.\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'ctrl-p'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : brings up the print dialog.\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'ctrl-f'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : brings up the \"Copy Frame Layout\" dialog.  This allows you to take the layout of\n"
			"\t\tone frame and apply it to one or more other frames.\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'ctrl-g'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : brings up the \"Copy Partition Layout\" dialog.  This allows you to take the layout of\n"
			"\t\tone partition and apply it to one or more other partitions.\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'ctrl-n'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : brings up a Dialog which allows you to create a new graph from a structure file.\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'ctrl-o'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : brings up a Dialog which allows you to open a previously saved graph (gvp file).\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'ctrl-s'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : saves the current graph.\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'ctrl-S'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : brings up the \"Save As\" dialog.\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'ctrl-w'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : closes the current graph.\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'ctrl-q'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : quits gmtkViz.\n");


	help_msg->SetDefaultStyle(bold);
	help_msg->AppendText("\nMouse Commands:\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'button 1'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : If the mouse pointer is on a selectable item (node, controlpoint, frame boarder, label) will select that item\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'button 1 + Shift'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : \n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'button 1 + drag motion'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : If drag is started on a selectable item then this will move the current selection\n"
		"\t\tIf drag is started elsewhere a box will be created (and destroyed when the mouse button is released)\n"
		"\t\teverything in this box will be selected\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'button 2'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : If the mouse pointer is on a control point will create a new control point\n"
		"\t\tIf there is only a strait line between nodes then a click on this line will create a spline and a control point\n");
	help_msg->SetDefaultStyle(italic);
	help_msg->AppendText("\t'button 3'");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText(" : If the mouse pointer is over a node will pop up node info\n");

	//help on the .Xdefaults File
	help_msg->SetDefaultStyle(bold);
	help_msg->AppendText("\n.Xdefaults File:\n");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText("\nThe .Xdefaults file allows you to set defaults for Pen Color, Width and Style as well Font Size and Style.");
	help_msg->AppendText("\nthe format is this:");
	help_msg->AppendText("\n\tgmtkViz.<PenName>Color: <ColorName>");
	help_msg->AppendText("\n\tgmtkViz.<PenName>Color: r,b,g (where r, g and b are integers between 0 and 255 inclusive)");
	help_msg->AppendText("\n\tgmtkViz.<PenName>Width: <integer>");
	help_msg->AppendText("\n\tgmtkViz.<PenName>Style: <StyleName>");

	help_msg->AppendText("\n\nPen Names are abreviated.\nValid PenNames are: ");

	//the interator for the PenDefaultMap
	map<const wxString, PenDefault_t *, ltstr>::iterator pen_iterator;
	//print the name of each pen
	unsigned int i = 0;
	for (pen_iterator = PenDefaultMap.begin(); pen_iterator != PenDefaultMap.end(); pen_iterator++){
		help_msg->AppendText(pen_iterator->first.c_str());
		if (i != PenDefaultMap.size() - 1){
			help_msg->AppendText(", ");
			if (i % 5 == 4)
				help_msg->AppendText("\n\t");
		}
		i++;
	}

	help_msg->AppendText("\n\nValid Colors are: "
			"\n\tAQUAMARINE, BLACK, BLUE, BLUE VIOLET, BROWN, CADET BLUE, CORAL,"
			"\n\tCORNFLOWER BLUE, CYAN, DARK GREY, DARK GREEN, DARK OLIVE GREEN, DARK ORCHID,"
			"\n\tDARK SLATE BLUE, DARK SLATE GREY DARK TURQUOISE, DIM GREY, FIREBRICK, FOREST"
			"\n\tGREEN, GOLD, GOLDENROD, GREY, GREEN, GREEN YELLOW, INDIAN RED, KHAKI, LIGHT"
			"\n\tBLUE, LIGHT GREY, LIGHT STEEL BLUE, LIME GREEN, MAGENTA, MAROON, MEDIUM"
			"\n\tAQUAMARINE, MEDIUM BLUE, MEDIUM FOREST GREEN, MEDIUM GOLDENROD, MEDIUM ORCHID,"
			"\n\tMEDIUM SEA GREEN, MEDIUM SLATE BLUE, MEDIUM SPRING GREEN, MEDIUM TURQUOISE,"
			"\n\tMEDIUM VIOLET RED, MIDNIGHT BLUE, NAVY, ORANGE, ORANGE RED, ORCHID, PALE GREEN,"
			"\n\tPINK, PLUM, PURPLE, RED, SALMON, SEA GREEN, SIENNA, SKY BLUE, SLATE BLUE,"
			"\n\tSPRING GREEN, STEEL BLUE, TAN, THISTLE, TURQUOISE, VIOLET, VIOLET RED, WHEAT,"
			"\n\tWHITE, YELLOW, YELLOW GREEN.");

	help_msg->AppendText("\n\nValid StyleName are: wxSOLID, wxDOT, wxLONG_DASH,");
   help_msg->AppendText(" wxSHORT_DASH, wxDOT_DASH, wxTRANSPARENT");

	help_msg->AppendText("\n\nFont Size is given in the format:");
	help_msg->AppendText("\n\tgmtkViz.fontSize: <num> (where num is an positive no zero integer)");

	help_msg->AppendText("\n\nFont Style defaults only work in version of wxWidgets 2.6.1 and beyond.");
	help_msg->AppendText("\nThe format for Font Style is:");
	help_msg->AppendText("\n\tgmtkViz.fontStyle: <style>");
	help_msg->AppendText("\nStyle is OS dependent, it is simply the name of a font, like Times New Roman, or Utopia.");
	help_msg->AppendText("\nTo find out what Styles are avalible for your system simply look at the font customization dialog in gmtkViz.");

	help_msg->AppendText("\n\nExamples:");
	help_msg->AppendText("\n\tgmtkViz.detPenWidth: 2");
	help_msg->AppendText("\n\tgmtkViz.detPenColor: 0,255,0");
	help_msg->AppendText("\n\tgmtkViz.randPenStyle: wxDOT");
	help_msg->AppendText("\n\tgmtkViz.randPenColor: VIOLET");
	help_msg->AppendText("\n\tgmtkViz.fontSize: 12");
	help_msg->AppendText("\n\tgmtkViz.fontName: Utopia");

	help_msg->SetDefaultStyle(bold);
	help_msg->AppendText("\n\nExtra Notes:\n");
	help_msg->SetDefaultStyle(normal);
	help_msg->AppendText("Once you make a selection invisible and you create a new selection or click on the canvas you can only make the\n"
		"nodes visible by toggling the visibilty of the nodes individually [by guessing their location] or by selecting\n"
		"View->\"Make All Nodes Visible\" from the tool bar\n"
		"\nIn the View Menu on the tool bar there are also options to Make All Node Labels Visible, Make All Frame Labels Visible\n"
		"and Hide the labels in the current selection\n");

	help_msg->AppendText("\nTo make the nodes in every frame except the frame under the mouse hit z x e\n"
		"this will select all visible nodes in the current frame, invert the selection and finally toggle the visibility of that selection\n");

	help_msg->AppendText("\n\nFor more detailed help go to: https://ssli.ee.washington.edu/ssliwiki/GMTK");

	help_sizer->Add(help_msg, 1, wxEXPAND, 0);

	help_sizer->Add( new wxButton( this, wxID_CLOSE, "Close" ), 0, wxALL | wxALIGN_CENTRE, 10 );

	SetAutoLayout(true);
	SetSizer(help_sizer);
	help_sizer->Fit(this);
	help_sizer->SetSizeHints(this);
	Layout();
}
