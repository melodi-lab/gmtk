// -*- C++ -*-

#include "GMTK_FileParser.h"
#include "GMTK_RVInfo.h"
#include <wx/wx.h>
#include <wx/dcbuffer.h>
#include <wx/ffile.h>
//#include <wx/gdicmn.h>
#include <wx/image.h>
#include <wx/list.h>
#include <wx/notebook.h>
#include <wx/print.h>
#include <wx/printdlg.h>
#include <wx/statline.h>
#include <wx/textfile.h>
#include <wx/tooltip.h>
#include <cassert>
#include <cmath>
#include <vector>

// Apprently these are needed in order to do anything with gmtk (even
// if you don't actually use them).
#include "GMTK_GMParms.h"
#include "rand.h"
#include "GMTK_ObservationMatrix.h"
GMParms GM_Parms;
RAND rnd(0);
ObservationMatrix globalObservationMatrix;

#define ACTUAL_SCALE (16)
static double gZoomMap[] = {
    0.0625/ACTUAL_SCALE,
    0.125/ACTUAL_SCALE,
    0.25/ACTUAL_SCALE,
    0.5/ACTUAL_SCALE,
    0.70710678118654746/ACTUAL_SCALE,
    0.84089641525371461/ACTUAL_SCALE,
    0.91700404320467122/ACTUAL_SCALE,
    1.0/ACTUAL_SCALE,
    1.0905077326652577/ACTUAL_SCALE,
    1.189207115002721/ACTUAL_SCALE,
    1.4142135623730951/ACTUAL_SCALE,
    2.0/ACTUAL_SCALE,
    4.0/ACTUAL_SCALE,
    8.0/ACTUAL_SCALE,
    16.0/ACTUAL_SCALE
};

#define NODE_RADIUS (10*ACTUAL_SCALE)
#define NEW_CP_OFFSET (10*ACTUAL_SCALE)
#define ARROW_LEN (12*ACTUAL_SCALE)
#define ARROW_WID (4*ACTUAL_SCALE)

class NameTag;
class VizNode;
class VizArc;
class VizSep;
class Selectable;
class ControlPoint;

class StructPage: public wxScrolledWindow
{
 public:
    StructPage() {}
    StructPage(wxWindow *parent, wxWindowID id,
	       wxFrame *parentFrame, wxNotebook *parentNotebook,
	       const wxString &file, bool old = true);
    ~StructPage();
    void OnPaint( wxPaintEvent &event );
    void OnChar( wxKeyEvent &event );
    void OnMouseEvent( wxMouseEvent &event );

    void Save( void );
    void SaveAs( void );
    bool RequestClose( void );

    bool Ready( void ) { return !gvpAborted; }

    void onComeForward( void );

    // things related to selecting
    Selectable *itemAt( const wxPoint& pt );
    void setAllSelected( bool newSelected );
    void setEndpointsSelected( bool newSelected, RVInfo::rvParent );
    void toggleSelectedInRect( const wxRect& rect );
    void moveSelected( int dx, int dy );

    int getWidth( void ) { return canvasWidth; }
    int getHeight( void ) { return canvasHeight; }
    int getScale( void ) { return displayScale; }
    void setScale( int newScale );
    void getName( wxString& name );
    void draw( wxDC& dc );
    wxPen switchingPen;
    wxPen conditionalPen;
    wxPen bothPen;
    wxPen frameBorderPen;
    wxPen chunkBorderPen;
    wxPen controlPointPen;
    wxPen nodePen;

    bool getViewCPs( void ) { return drawCPs; }
    bool getViewLines( void ) { return drawLines; }
    bool getViewSplines( void ) { return drawSplines; }
    bool getViewArrowHeads( void ) { return drawArrowHeads; }
    bool getViewNodes( void ) { return drawNodes; }
    bool getViewDirectLines( void ) { return drawDirectLines; }
    bool getViewFrameSeps( void ) { return drawFrameSeps; }
    bool getViewNodeNames( void ) { return drawNodeNames; }
    bool getViewToolTips( void ) { return drawToolTips; }
    void toggleViewCPs( void );
    void toggleViewLines( void );
    void toggleViewSplines( void );
    void toggleViewArrowHeads( void );
    void toggleViewNodes( void );
    void toggleViewDirectLines( void );
    void toggleViewFrameSeps( void );
    void toggleViewNodeNames( void );
    void toggleViewToolTips( void );

    DECLARE_DYNAMIC_CLASS(StructPage)
    DECLARE_EVENT_TABLE()
private:
    wxFrame *parentFrame;
    wxNotebook *parentNotebook;
    int pageNum;
    wxBitmap *content; // the drawing buffer
    std::map< RVInfo::rvParent, unsigned int > nameVizNodeMap;
    wxRect selectBox;
    bool drawSelectBox;
    bool drawCPs;
    bool drawLines;
    bool drawSplines;
    bool drawArrowHeads;
    bool drawNodes;
    bool drawDirectLines;
    bool drawFrameSeps;
    bool drawNodeNames;
    bool drawToolTips;
    int displayScale;
    long canvasWidth;
    long canvasHeight;
    void initNodes( void );
    void initArcs( void );
    void redraw( void );
    void blit( wxDC& dc );
    void blit( void );
    VizArc* newArc(int i, int j);
    VizArc* findArcOwning( ControlPoint *cp, int& index );
    void deleteSelectedCps( void );
    bool crossesNode( wxCoord x0, wxCoord x1 );
    bool crossesFrameEnd( wxCoord x0, wxCoord x1 );
    bool inBounds( wxCoord x, wxCoord y );

    // things related to the structure (str) file
    wxString strFile;
    FileParser *fp;

    // things related to the position (gvp) file
    wxString gvpFile;
    bool gvpDirty;
    bool gvpAborted;
    void fillMap( void );
    struct ltstr {
	bool operator()(const wxString s1, const wxString s2) const
	{
	    return (s1 < s2);
	}
    };
    map<const wxString, wxString, ltstr> config;
    std::vector< VizNode* > nodes; // node positions
    std::vector< NameTag* > nodeNameTags; // node nametags
    std::vector< std::vector< VizArc* > > arcs; // arcs[a][b] is the
    // information about the arc/spline from node a to node b (or NULL)
    std::vector< VizSep* > frameEnds; // positions of dividers between frames
};


// Represents anything that can be selected
class Selectable {
public:
    virtual bool getSelected( void ) { return selected; }
    virtual void setSelected( bool newSelected ) { selected = newSelected; }
    virtual void toggleSelected( void ) { setSelected(!getSelected()); }
    virtual bool onMe( const wxPoint& pt ) = 0;
    virtual bool inRect( const wxRect& rect ) = 0;
    		Selectable() { selected = false; }
    virtual 	~Selectable() {  }
protected:
    bool selected;
};


// Represents a node's nametag
class NameTag : public Selectable {
public:
    wxPoint pos;
    wxPoint size;
    wxString name;
    NameTag( const wxPoint& newPos, const wxString& newName );
    void draw( wxDC *dc );
    virtual bool onMe( const wxPoint& pt );
    virtual bool inRect( const wxRect& rect );
};

// Represents a node 
class VizNode : public Selectable {
public:
    wxPoint center;
    RVInfo *rvi;
    StructPage *page;
    RVInfo::rvParent rvId;
    NameTag nametag;
    wxWindow *tipWin;
    VizNode( const wxPoint& newPos, RVInfo *newRvi, StructPage *parentPage );
    ~VizNode( void );
    void draw( wxDC *dc );
    // Selectable methods
    virtual void setSelected( bool newSelected );
    virtual bool onMe( const wxPoint& pt );
    virtual bool inRect( const wxRect& rect ) { return rect.Inside(center); }
};


// Represents a control point
class ControlPoint : public Selectable {
public:
    wxPoint pos;
    VizArc *arc;
    ControlPoint( const wxPoint& pt );
    virtual bool onMe( const wxPoint& pt );
    virtual bool inRect( const wxRect& rect ) { return rect.Inside(pos); }
};

// Represents an arc
class VizArc {
public:
    wxList *points;
    std::vector< ControlPoint* > *cps;
    StructPage *page;
    bool switching;
    bool conditional;
    VizArc( std::vector< ControlPoint* > *newCps, StructPage *newPage );
    ~VizArc( void );
    void draw( wxDC *dc );
};


// Represents a separator between frames
class VizSep : public Selectable {
public:
    wxCoord x;
    bool chunkBorder;
    StructPage *page;
    VizSep( wxCoord xNew, StructPage *newPage, bool newChunkBorder = false );
    void draw( wxDC *dc );
    virtual bool onMe( const wxPoint& pt );
    virtual bool inRect( const wxRect& rect ) { return false; }
};


class GmtkPrintout : public wxPrintout {
public:
    GmtkPrintout(StructPage *newPage, const wxChar *title = wxT("gmtkViz printout"));
    bool OnPrintPage(int page);
    bool HasPage(int page);
    void GetPageInfo(int *minPage, int *maxPage, int *selPageFrom, int *selPageTo);
    void DrawPageOne(wxDC *dc);
private:
    StructPage *page;
};

class GFrame: public wxFrame {
public:
    enum {
        MENU_FILE_NEW = 1006,
        MENU_FILE_OPEN,
        MENU_FILE_SAVE,
        MENU_FILE_SAVEAS,
        MENU_FILE_PAGESETUP,
        MENU_FILE_PRINT,
        MENU_FILE_CLOSE,
        MENU_FILE_EXIT,
	MENU_VIEW_CPS,
	MENU_VIEW_LINES,
	MENU_VIEW_SPLINES,
	MENU_VIEW_ARROW_HEADS,
	MENU_VIEW_NODES,
	MENU_VIEW_DIRECT_LINES,
	MENU_VIEW_FRAME_SEPS,
	MENU_VIEW_NODE_NAMES,
	MENU_VIEW_TOOLTIPS,
	MENU_ZOOM_BEGIN,
	MENU_ZOOM_2_pow_neg_4dot000,
	MENU_ZOOM_2_pow_neg_3dot000,
	MENU_ZOOM_2_pow_neg_2dot000,
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
	MENU_ZOOM_2_pow_pos_3dot000,
	MENU_ZOOM_2_pow_pos_4dot000,
	MENU_ZOOM_END
    };

    /**
     * Constructor. Creates a new GFrame.
     */
    GFrame(wxWindow* parent, int id, const wxString& title, const wxPoint& pos=wxDefaultPosition, const wxSize& size=wxDefaultSize, long style=wxDEFAULT_FRAME_STYLE);

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

    void OnMenuViewCPs(wxCommandEvent &event);
    void OnMenuViewLines(wxCommandEvent &event);
    void OnMenuViewSplines(wxCommandEvent &event);
    void OnMenuViewArrowHeads(wxCommandEvent &event);
    void OnMenuViewNodes(wxCommandEvent &event);
    void OnMenuViewDirectLines(wxCommandEvent &event);
    void OnMenuViewFrameSeps(wxCommandEvent &event);
    void OnMenuViewNodeNames(wxCommandEvent &event);
    void OnMenuViewToolTips(wxCommandEvent &event);
    void OnMenuZoom(wxCommandEvent &event);
    void OnNotebookPageChanged(wxCommandEvent &event);

private:
    void set_properties();
    void do_layout();
    wxPrintData printData;
    wxPageSetupData pageSetupData;

protected:
    wxMenuBar* MainVizWindow_menubar;
    wxStatusBar* MainVizWindow_statusbar;
    wxStaticText* about_label;
    wxTextCtrl* about_info;
    wxPanel* about_pane;
    wxNotebook* struct_notebook;
    DECLARE_EVENT_TABLE()
};

class GMTKStructVizApp: public wxApp {
public:
    bool OnInit();
};

/*** implementation starts here ***/

IMPLEMENT_APP(GMTKStructVizApp)

bool GMTKStructVizApp::OnInit()
{
    wxInitAllImageHandlers();
    GFrame* MainVizWindow = new GFrame(0, -1, "");
    SetTopWindow(MainVizWindow);
    MainVizWindow->Show();
    return true;
}

GFrame::GFrame( wxWindow* parent, int id, const wxString& title,
		const wxPoint& pos, const wxSize& size, long style )
    : wxFrame( parent, id, title, pos, size, wxDEFAULT_FRAME_STYLE )
{
    struct_notebook = new wxNotebook(this, -1, wxDefaultPosition, wxDefaultSize, wxCLIP_CHILDREN);
    about_pane = new wxPanel(struct_notebook, -1);
    MainVizWindow_menubar = new wxMenuBar();
    SetMenuBar(MainVizWindow_menubar);
    wxMenu* menu_file = new wxMenu();
    menu_file->Append(MENU_FILE_NEW, wxT("&New...\tCtrl+N"), wxT("Create a new placement file (requires an existing structure file)"), wxITEM_NORMAL);
    menu_file->Append(MENU_FILE_OPEN, wxT("&Open...\tCtrl+O"), wxT("Open an existing placement file"), wxITEM_NORMAL);
    menu_file->AppendSeparator();
    menu_file->Append(MENU_FILE_SAVE, wxT("&Save\tCtrl+S"), wxT("Save the current placement file"), wxITEM_NORMAL);
    menu_file->Append(MENU_FILE_SAVEAS, wxT("Save &As...\tCtrl+Shift+S"), wxT("Save the current placement file with a different name"), wxITEM_NORMAL);
    menu_file->AppendSeparator();
    menu_file->Append(MENU_FILE_PAGESETUP, wxT("Page Setup..."), wxT("Set up page size/orientation"), wxITEM_NORMAL);
    menu_file->Append(MENU_FILE_PRINT, wxT("&Print...\tCtrl+P"), wxT("Preview and print the current graph"), wxITEM_NORMAL);
    menu_file->AppendSeparator();
    menu_file->Append(MENU_FILE_CLOSE, wxT("&Close\tCtrl+W"), wxT("Close current placement file"), wxITEM_NORMAL);
    menu_file->Append(MENU_FILE_EXIT, wxT("E&xit\tCtrl+Q"), wxT("Close all files and exit"), wxITEM_NORMAL);
    MainVizWindow_menubar->Append(menu_file, wxT("&File"));
    wxMenu* menu_view = new wxMenu();
    menu_view->Append(MENU_VIEW_CPS, wxT("Draw Control Points"), wxT("Toggle display of arc spline control points"), wxITEM_CHECK);
    menu_view->Append(MENU_VIEW_LINES, wxT("Draw Arc Lines"), wxT("Toggle display of straight lines between control points in arcs"), wxITEM_CHECK);
    menu_view->Append(MENU_VIEW_SPLINES, wxT("Draw Arc Splines"), wxT("Toggle display of arc splines"), wxITEM_CHECK);
    menu_view->Append(MENU_VIEW_ARROW_HEADS, wxT("Draw Arrow Heads"), wxT("Toggle display of arrow heads on arcs"), wxITEM_CHECK);
    menu_view->Append(MENU_VIEW_NODES, wxT("Draw Nodes"), wxT("Toggle display of nodes"), wxITEM_CHECK);
    menu_view->Append(MENU_VIEW_DIRECT_LINES, wxT("Draw Direct Lines"), wxT("Toggle display of direct straight lines for arcs"), wxITEM_CHECK);
    menu_view->Append(MENU_VIEW_FRAME_SEPS, wxT("Draw Frame Separators"), wxT("Toggle display of frame separators"), wxITEM_CHECK);
    menu_view->Append(MENU_VIEW_NODE_NAMES, wxT("Draw Node Names"), wxT("Toggle display of node names"), wxITEM_CHECK);
    // XXX: menu_view->Append(MENU_VIEW_TOOLTIPS, wxT("Draw Tool Tips"), wxT("Toggle display of tool tips for node names"), wxITEM_CHECK);
    MainVizWindow_menubar->Append(menu_view, wxT("View"));
    MainVizWindow_menubar->Enable(MENU_FILE_SAVE, false);
    MainVizWindow_menubar->Enable(MENU_FILE_SAVEAS, false);
    MainVizWindow_menubar->Enable(MENU_FILE_PRINT, false);
    MainVizWindow_menubar->Enable(MENU_FILE_CLOSE, false);
    MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("View")), false);
    wxMenu* menu_zoom = new wxMenu();
    for (int i = 0; i < MENU_ZOOM_END - MENU_ZOOM_BEGIN - 1; i++) {
	wxString zoomStr;
	zoomStr.sprintf("%7.3f", gZoomMap[i]*(ACTUAL_SCALE));
	menu_zoom->Append( i + MENU_ZOOM_BEGIN + 1, zoomStr,
			   wxEmptyString, wxITEM_RADIO );
    }
    MainVizWindow_menubar->Append(menu_zoom, wxT("Zoom"));
    MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("Zoom")), false);
    MainVizWindow_statusbar = CreateStatusBar(2);
    about_label = new wxStaticText(about_pane, -1, wxT("GMTKStructViz version 0.0.1\n\nThe graph visualizer and organizer for GMTK structure files.\n\nwritten by Evan Dower <evantd@ssli.ee.washington.edu>"), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE);
    about_info = new wxTextCtrl(about_pane, -1, wxT("GMTKStructViz reads GMTK structure files and attempts to display their contents semi-intelligently. Since it is not human, it can only have limited success in this domain. Thus the user is permitted to move nodes and edges to organize the graph in a more logical and visually appealing way than GMTKStructViz's original guess."), wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE|wxTE_READONLY);

    set_properties();
    do_layout();
}


void GFrame::set_properties()
{
    SetTitle(wxT("GMTK Structure File Visualizer"));
    int MainVizWindow_statusbar_widths[] = { -1, 16 };
    MainVizWindow_statusbar->SetStatusWidths(2, MainVizWindow_statusbar_widths);
    const wxString MainVizWindow_statusbar_fields[] = {
	wxT("About GMTKStructViz"),
        wxEmptyString
    };
    for(int i = 0; i < MainVizWindow_statusbar->GetFieldsCount(); ++i) {
        MainVizWindow_statusbar->SetStatusText(MainVizWindow_statusbar_fields[i], i);
    }
}


void GFrame::do_layout()
{
    wxBoxSizer* MainVizWindow_sizer = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer* about_sizer = new wxBoxSizer(wxVERTICAL);
    wxStaticBoxSizer* about_label_static_sizer = new wxStaticBoxSizer(new wxStaticBox(about_pane, -1, wxT("")), wxHORIZONTAL);
    about_label_static_sizer->Add(about_label, 0, wxALL|wxALIGN_CENTER_HORIZONTAL, 2);
    about_sizer->Add(about_label_static_sizer, 1, wxEXPAND, 0);
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

BEGIN_EVENT_TABLE(GFrame, wxFrame)
    EVT_MENU(MENU_FILE_NEW, GFrame::OnMenuFileNew)
    EVT_MENU(MENU_FILE_OPEN, GFrame::OnMenuFileOpen)
    EVT_MENU(MENU_FILE_SAVE, GFrame::OnMenuFileSave)
    EVT_MENU(MENU_FILE_SAVEAS, GFrame::OnMenuFileSaveas)
    EVT_MENU(MENU_FILE_PAGESETUP, GFrame::OnMenuFilePageSetup)
    EVT_MENU(MENU_FILE_PRINT, GFrame::OnMenuFilePrint)
    EVT_MENU(MENU_FILE_CLOSE, GFrame::OnMenuFileClose)
    EVT_MENU(MENU_FILE_EXIT, GFrame::OnMenuFileExit)
    EVT_MENU(MENU_VIEW_CPS, GFrame::OnMenuViewCPs)
    EVT_MENU(MENU_VIEW_LINES, GFrame::OnMenuViewLines)
    EVT_MENU(MENU_VIEW_SPLINES, GFrame::OnMenuViewSplines)
    EVT_MENU(MENU_VIEW_ARROW_HEADS, GFrame::OnMenuViewArrowHeads)
    EVT_MENU(MENU_VIEW_NODES, GFrame::OnMenuViewNodes)
    EVT_MENU(MENU_VIEW_DIRECT_LINES, GFrame::OnMenuViewDirectLines)
    EVT_MENU(MENU_VIEW_FRAME_SEPS, GFrame::OnMenuViewFrameSeps)
    EVT_MENU(MENU_VIEW_NODE_NAMES, GFrame::OnMenuViewNodeNames)
    EVT_MENU(MENU_VIEW_TOOLTIPS, GFrame::OnMenuViewToolTips)
    EVT_MENU_RANGE(MENU_ZOOM_BEGIN+1, MENU_ZOOM_END-1, GFrame::OnMenuZoom)
    EVT_NOTEBOOK_PAGE_CHANGED(wxID_ANY, GFrame::OnNotebookPageChanged)
    EVT_CLOSE(GFrame::OnClose)
END_EVENT_TABLE()

void GFrame::OnMenuFileNew(wxCommandEvent &event)
{
    wxFileDialog dlg(this,
		     "Find the desired structure file", "", "",
		     "All files|*"
		     "|Structure Files|*.str"
		     "|Text Files|*.txt;*.text",
		     wxOPEN | wxCHANGE_DIR, wxDefaultPosition);

    dlg.SetFilterIndex(1); // show only .str's by default

    if ( dlg.ShowModal() == wxID_OK ) {
	new StructPage(struct_notebook, -1, this, struct_notebook,
		       dlg.GetPath(), false);
	wxCommandEvent dummy;
	OnNotebookPageChanged(dummy);
    }
}

void GFrame::OnMenuFileOpen(wxCommandEvent &event)
{
    wxFileDialog dlg(this,
		     "Find the desired position file", "", "",
		     "All files|*"
		     "|Position Files|*.gvp"
		     "|Text Files|*.txt;*.text",
		     wxOPEN | wxCHANGE_DIR, wxDefaultPosition);

    dlg.SetFilterIndex(1); // show only .gvp's by default

    if ( dlg.ShowModal() == wxID_OK ) {
	/* demonstrate that we got the right file
	wxTextCtrl *new_text_area = new wxTextCtrl(struct_notebook, -1,
						   wxEmptyString,
						   wxDefaultPosition,
						   wxDefaultSize,
						   wxTE_MULTILINE);
	new_text_area->LoadFile(dlg.GetPath());
	struct_notebook->AddPage(new_text_area, dlg.GetFilename());
	SetStatusText(dlg.GetFilename(), 0);*/

	// This will add itself to the notebook and be destroyed when
	// the notebook is destroyed.
	new StructPage( struct_notebook, -1, this, struct_notebook,
			dlg.GetPath(), true );
	wxCommandEvent dummy;
	OnNotebookPageChanged(dummy);
    }
}

void GFrame::OnMenuFileSave(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage)
	curPage->Save();
}

void GFrame::OnMenuFileSaveas(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage)
	curPage->SaveAs();
}

void GFrame::OnMenuFilePageSetup(wxCommandEvent &event)
{
    pageSetupData = printData;

    wxPageSetupDialog pageSetupDialog(this, &pageSetupData);

    pageSetupDialog.ShowModal();

    printData = pageSetupDialog.GetPageSetupData().GetPrintData();
    pageSetupData = pageSetupDialog.GetPageSetupData();
}

void GFrame::OnMenuFilePrint(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage) {
	wxString name;
	curPage->getName(name);
	// Pass two printout objects: for preview, and possible printing.
	wxPrintDialogData printDialogData(printData);
	wxPrintPreview *preview =
	    new wxPrintPreview(new GmtkPrintout(curPage, name),
			       new GmtkPrintout(curPage, name),
			       &printDialogData);
	if (!preview->Ok()) {
	    delete preview;
	    wxMessageBox(_T("There was a problem previewing.\nPerhaps your current printer is not set correctly?"), wxT("Previewing"), wxOK);
	    return;
	}
	
	wxPreviewFrame *frame =
	    new wxPreviewFrame(preview, this, wxT("gmtkViz Print Preview"),
			       wxPoint(100, 100), wxSize(600, 650));
	frame->Centre(wxBOTH);
	frame->Initialize();
	frame->Show(TRUE);
    }
}

void GFrame::OnMenuFileClose(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage) {
	if (curPage->RequestClose()) {
	    struct_notebook->DeletePage(curPageNum);
	}
    }
}

void GFrame::OnMenuFileExit(wxCommandEvent &event)
{
    Close(false);
}

void
GFrame::OnClose(wxCloseEvent& event)
{
    bool destroy = true;

    if ( event.CanVeto() ) {
	StructPage *curPage;
	for ( int pageNum = struct_notebook->GetPageCount() - 1;
	      pageNum >= 0 && destroy;
	      pageNum-- ) {
	    curPage = dynamic_cast<StructPage*>
		(struct_notebook->GetPage(pageNum));
	    if (curPage) {
		if (curPage->RequestClose()) {
		    struct_notebook->DeletePage(pageNum);
		} else {
		    event.Veto();
		    destroy = false;
		}
	    }
	}
    }
    if ( destroy ) {
	Destroy();
    }
}

void
GFrame::OnMenuViewCPs(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage)
	curPage->toggleViewCPs();
}

void
GFrame::OnMenuViewLines(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage)
	curPage->toggleViewLines();
}

void
GFrame::OnMenuViewSplines(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage)
	curPage->toggleViewSplines();
}

void
GFrame::OnMenuViewArrowHeads(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage)
	curPage->toggleViewArrowHeads();
}

void
GFrame::OnMenuViewNodes(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage)
	curPage->toggleViewNodes();
}

void
GFrame::OnMenuViewDirectLines(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage)
	curPage->toggleViewDirectLines();
}

void
GFrame::OnMenuViewFrameSeps(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage)
	curPage->toggleViewFrameSeps();
}

void
GFrame::OnMenuViewNodeNames(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage)
	curPage->toggleViewNodeNames();
}

void
GFrame::OnMenuViewToolTips(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage)
	curPage->toggleViewToolTips();
}

void
GFrame::OnMenuZoom(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage) {
	int id = event.GetId()/*, scale = curPage->getScale()*/;
	//MainVizWindow_menubar->Check(MENU_ZOOM_BEGIN + scale + 1, false);
	MainVizWindow_menubar->Check(id, true);
	curPage->setScale(id - MENU_ZOOM_BEGIN - 1);
    }
}

void
GFrame::OnNotebookPageChanged(wxCommandEvent &event)
{
    // figure out which page this is for and pass the buck
    int curPageNum = struct_notebook->GetSelection();
    StructPage *curPage = dynamic_cast<StructPage*>
	(struct_notebook->GetPage(curPageNum));
    if (curPage) {
	curPage->onComeForward();
	MainVizWindow_menubar->Enable(MENU_FILE_SAVE, true);
	MainVizWindow_menubar->Enable(MENU_FILE_SAVEAS, true);
	MainVizWindow_menubar->Enable(MENU_FILE_PRINT, true);
	MainVizWindow_menubar->Enable(MENU_FILE_CLOSE, true);
	MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("View")), true);
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
	MainVizWindow_menubar->Check( MENU_VIEW_DIRECT_LINES,
				      curPage->getViewDirectLines() );
	MainVizWindow_menubar->Check( MENU_VIEW_FRAME_SEPS,
				      curPage->getViewFrameSeps() );
	MainVizWindow_menubar->Check( MENU_VIEW_NODE_NAMES,
				      curPage->getViewNodeNames());
	// XXX: MainVizWindow_menubar->Check( MENU_VIEW_TOOLTIPS, curPage->getViewToolTips() );
	MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("Zoom")), true);
	/*for (int i = 0, scale = curPage->getScale(); i < MENU_ZOOM_END - MENU_ZOOM_BEGIN - 1; i++) {
	    MainVizWindow_menubar->Check( i + MENU_ZOOM_BEGIN + 1,
					  i==scale );
	}*/
	MainVizWindow_menubar->Check( curPage->getScale()+MENU_ZOOM_BEGIN+1,
				      true );
    }
    else {
	SetStatusText(wxEmptyString, 1);
	SetStatusText(wxT("About GMTKStructViz"), 0);
	MainVizWindow_menubar->Enable(MENU_FILE_SAVE, false);
	MainVizWindow_menubar->Enable(MENU_FILE_SAVEAS, false);
	MainVizWindow_menubar->Enable(MENU_FILE_PRINT, false);
	MainVizWindow_menubar->Enable(MENU_FILE_CLOSE, false);
	MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("View")), false);
	MainVizWindow_menubar->EnableTop(MainVizWindow_menubar->FindMenu(wxT("Zoom")), false);
    }
}

IMPLEMENT_DYNAMIC_CLASS(StructPage, wxScrolledWindow)

BEGIN_EVENT_TABLE(StructPage, wxScrolledWindow)
    EVT_PAINT(StructPage::OnPaint)
    EVT_MOUSE_EVENTS(StructPage::OnMouseEvent)
    EVT_CHAR(StructPage::OnChar)
END_EVENT_TABLE()

StructPage::StructPage(wxWindow *parent, wxWindowID id,
		       wxFrame *parentFrame, wxNotebook *parentNotebook,
		       const wxString &file, bool old)
    : wxScrolledWindow( parent, id, wxDefaultPosition, wxDefaultSize,
			wxSUNKEN_BORDER | wxTAB_TRAVERSAL, _T("") ),
      switchingPen(*wxCYAN_PEN), conditionalPen(*wxBLACK_PEN),
      bothPen(*wxRED_PEN), frameBorderPen(*wxLIGHT_GREY_PEN),
      chunkBorderPen(*wxBLACK_PEN), controlPointPen(*wxRED_PEN),
      nodePen(*wxBLACK_PEN)
{
    this->parentFrame = parentFrame;
    this->parentNotebook = parentNotebook;
    gvpDirty = false;
    gvpAborted = false;
    drawSelectBox = false;
    drawCPs = true;
    drawLines = false;
    drawSplines = true;
    drawArrowHeads = true;
    drawNodes = true;
    drawDirectLines = false;
    drawFrameSeps = true;
    drawNodeNames = true;
    drawToolTips = true;
    content = NULL;

    //switchingPen.SetStyle(wxDOT);
    //conditionalPen.SetStyle(wxSOLID);
    //bothPen.SetStyle(wxLONG_DASH);
    //frameBorderPen.SetStyle(wxDOT);
    //chunkBorderPen.SetStyle(wxSOLID);
    switchingPen.SetWidth(ACTUAL_SCALE);
    conditionalPen.SetWidth(ACTUAL_SCALE);
    bothPen.SetWidth(ACTUAL_SCALE);
    frameBorderPen.SetWidth(ACTUAL_SCALE);
    chunkBorderPen.SetWidth(ACTUAL_SCALE);
    controlPointPen.SetWidth(ACTUAL_SCALE);
    controlPointPen.SetJoin(wxJOIN_MITER);
    nodePen.SetWidth(ACTUAL_SCALE);
    SetScrollRate( 10, 10 );
    wxToolTip::SetDelay(250);
    wxToolTip::Enable(drawToolTips);

    // mangle the name and set the tab title and status bar
    int beginning, ending;
    if ( (beginning = file.rfind('/') + 1) < 1 )
	beginning = 0;
    if ( (ending = file.rfind(_T(".str"))) <= beginning &&
	 (ending = file.rfind(_T(".gvp"))) <= beginning )
	ending = file.size();
    parentNotebook->InsertPage(pageNum = parentNotebook->GetPageCount(),
			       this, file.substr(beginning, ending), true);

    if (old) {
	gvpFile = file;
	fillMap();
    } else {
	strFile = file;
    }

    // load up the structure file
    fp = new FileParser(strFile, NULL);

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
    initArcs();
    // At this point we no longer need the file parser, and since it
    // uses global variables, we can't have more than one. Thus, we
    // delete it now, rather than later.
    delete fp;
    fp = NULL;
    if (gvpAborted)
	wxLogMessage("Position file contained errors: proceed with caution");

    onComeForward();

    setScale(7);
}

StructPage::~StructPage( void )
{
    int numNodes = nodes.size();

    for (int i = frameEnds.size() - 1; i >= 0; i--) {
	delete frameEnds[i];
	frameEnds[i] = NULL;
    }
    for (int i = numNodes - 1; i >= 0; i--) {
	for (int j = numNodes - 1; j >= 0; j--) {
	    if (arcs[i][j]) {
		delete arcs[i][j];
		arcs[i][j] = NULL;
	    }
	}
    }
    for (int i = numNodes - 1; i >= 0; i--) {
	delete nodes[i];
	nodes[i] = NULL;
    }
    delete content;
    content = NULL;
}

void
StructPage::fillMap( void )
{
    wxTextFile gvp(gvpFile);
    wxString line;

    if (gvp.Open()) {
	for (line = gvp.GetFirstLine(); !gvp.Eof(); line = gvp.GetNextLine()) {
	    config[line.BeforeFirst('=')] = line.AfterFirst('=');
	}
	strFile = config[wxT("strFile")];
	if (!strFile.length()) {
	    wxLogMessage(wxT("gvp file doesn't specify a structure file"));
	    gvpAborted = true;
	}
    } else {
	wxLogMessage(wxT("Failed to open position file for reading"));
	gvpAborted = true;
    }
}

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
    int numRows = (int)ceil(sqrt(2.0*numVars));
    int yMax = 0;
    wxString key, value;

    key.sprintf("numNodes");
    value = config[key];
    if (value != wxEmptyString) {
	long numNodes;
	if (value.ToLong(&numNodes)) {
	    if (numNodes != numVars) {
		wxLogMessage("gvp and str disagree on number of nodes");
		gvpAborted = true;
	    }
	} else {
	    wxLogMessage("numNodes is not a number");
	    gvpAborted = true;
	}
    }
    key.sprintf("numFrames");
    value = config[key];
    if (value != wxEmptyString) {
	long numFrames;
	if (value.ToLong(&numFrames)) {
	    if (numFrames != (long)fp->_maxFrame + 1) {
		wxLogMessage("gvp and str disagree on number of frames");
		gvpAborted = true;
	    }
	} else {
	    wxLogMessage("numFrames is not a number");
	    gvpAborted = true;
	}
    }

    for (int i = 0, row = 0; i < numVars; i++, row++) {
	if (row >= numRows || fp->rvInfoVector[i].frame != curFrame) {
	    if (curPos.y > yMax)
		yMax = curPos.y;
	    curPos.x += 80*ACTUAL_SCALE;
	    curPos.y = 120*ACTUAL_SCALE;
	    row = 0;
	}
	if (fp->rvInfoVector[i].frame != curFrame) {
	    frameEnds.push_back(new VizSep( curPos.x, this,
					    curFrame==fp->_firstChunkframe-1 ||
					    curFrame==fp->_lastChunkframe ));
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
	    curPos.x += 80*ACTUAL_SCALE;
	    curFrame++;
	    assert(fp->rvInfoVector[i].frame == curFrame);
	}
	curPos.y += 80*ACTUAL_SCALE;
	nodes.push_back(new VizNode( curPos, &fp->rvInfoVector[i], this ));
	assert(nodes.size() == (unsigned)i + 1);
	nodeNameTags.push_back(&nodes[i]->nametag);
	assert(nodeNameTags.size() == nodes.size());
	nameVizNodeMap[nodes[i]->rvId] = i;

	// if a position was specified, move the node to it
	// x
	key.sprintf("nodes[%d].center.x", i);
	value = config[key];
	if (value != wxEmptyString) {
	    long xPos;
	    if (value.ToLong(&xPos)) {
		nodes[i]->center.x = xPos;
		nodeNameTags[i]->pos.x = xPos;
	    } else {
		wxString msg;
		msg.sprintf("nodes[%d].center.x is not a number", i);
		wxLogMessage(msg);
		gvpAborted = true;
	    }
	}
	// y
	key.sprintf("nodes[%d].center.y", i);
	value = config[key];
	if (value != wxEmptyString) {
	    long yPos;
	    if (value.ToLong(&yPos)) {
		nodes[i]->center.y = yPos;
		nodeNameTags[i]->pos.y = yPos;
	    } else {
		wxString msg;
		msg.sprintf("nodes[%d].center.y is not a number", i);
		wxLogMessage(msg);
		gvpAborted = true;
	    }
	}
	// and move the nametag, if requested
	// x
	key.sprintf("nodes[%d].nametag.pos.x", i);
	value = config[key];
	if (value != wxEmptyString) {
	    long xPos;
	    if (value.ToLong(&xPos)) {
		nodeNameTags[i]->pos.x = xPos;
	    } else {
		wxString msg;
		msg.sprintf("nodes[%d].nametag.pos.x is not a number", i);
		wxLogMessage(msg);
		gvpAborted = true;
	    }
	}
	// y
	key.sprintf("nodes[%d].nametag.pos.y", i);
	value = config[key];
	if (value != wxEmptyString) {
	    long yPos;
	    if (value.ToLong(&yPos)) {
		nodeNameTags[i]->pos.y = yPos;
	    } else {
		wxString msg;
		msg.sprintf("nodes[%d].nametag.pos.y is not a number", i);
		wxLogMessage(msg);
		gvpAborted = true;
	    }
	}
    }
    if (curPos.y > yMax)
	yMax = curPos.y;
    canvasWidth = curPos.x + 200*ACTUAL_SCALE;
    canvasHeight = yMax + 200*ACTUAL_SCALE;
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
    gvpDirty = true;
    onComeForward();
}

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

    for (unsigned int j = 0; j < numNodes; j++) {
	// for each switching parent
	for ( unsigned int k = 0;
	      k < nodes[j]->rvi->switchingParents.size();
	      k++ ) {
	    RVInfo::rvParent absId;
	    absId.first = nodes[j]->rvi->switchingParents[k].first;
	    absId.second = nodes[j]->rvi->switchingParents[k].second
		+ nodes[j]->rvi->frame;
	    int i = nameVizNodeMap[ absId ];
	    if (!arcs[i][j])
		arcs[i][j] = newArc(i, j);
	    /*if (arcs[i][j]->switching) {
		wxString msg;
		msg.sprintf("switching arc from node %d (%s:%d) "
			    "to node %d (%s:%d) already exists",
			    i, nodes[i]->rvId.first.c_str(),
			    nodes[i]->rvId.second,
			    j, nodes[j]->rvId.first.c_str(),
			    nodes[j]->rvId.second);
		wxLogMessage(msg);
	    }*/
	    arcs[i][j]->switching = true;
	}
	// for each conditional parent
	for ( unsigned int k = 0;
	      k < nodes[j]->rvi->conditionalParents.size();
	      k++ ) {
	    for ( unsigned int l = 0;
		  l < nodes[j]->rvi->conditionalParents[k].size();
		  l++ ) {
		RVInfo::rvParent absId;
		absId.first = nodes[j]->rvi->conditionalParents[k][l].first;
		absId.second = nodes[j]->rvi->conditionalParents[k][l].second
		    + nodes[j]->rvi->frame;
		int i=nameVizNodeMap[absId];
		if (!arcs[i][j])
		    arcs[i][j] = newArc(i, j);
		/*if (arcs[i][j]->conditional) {
		    wxString msg;
		    msg.sprintf("conditional arc from node %d (%s:%d) "
				"to node %d (%s:%d) already exists",
				i, nodes[i]->rvId.first.c_str(),
				nodes[i]->rvId.second,
				j, nodes[j]->rvId.first.c_str(),
				nodes[j]->rvId.second);
		    wxLogMessage(msg);
		}*/
		arcs[i][j]->conditional = true;
	    }
	}
    }
    gvpDirty = true;
    onComeForward();
}

VizArc*
StructPage::newArc( int i, int j )
{
    std::vector< ControlPoint* > *cps
	= new std::vector< ControlPoint* >;
    wxString key, value;

    cps->clear();
    cps->push_back(new ControlPoint(nodes[i]->center));

    wxPoint pt;
    key.sprintf("arcs[%d][%d].numCPs", i, j);
    value = config[key];
    if (value != wxEmptyString) {
	long numCPs;
	if (value.ToLong(&numCPs)) {
	    for (int k = 1; k < numCPs - 1; k++) {
		key.sprintf("arcs[%d][%d].cps[%d].pos.x", i, j, k);
		value = config[key];
		if (value != wxEmptyString) {
		    long xPos;
		    if (value.ToLong(&xPos)) {
			key.sprintf("arcs[%d][%d].cps[%d].pos.y", i, j, k);
			value = config[key];
			if (value != wxEmptyString) {
			    long yPos;
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
	pt.x = (2*nodes[i]->center.x + nodes[j]->center.x)/3 + 20*ACTUAL_SCALE;
	pt.y = (2*nodes[i]->center.y + nodes[j]->center.y)/3 - 20*ACTUAL_SCALE;
	ControlPoint *mid = new ControlPoint(pt);
	cps->push_back(mid);
	pt.x = (nodes[i]->center.x + 2*nodes[j]->center.x)/3 + 20*ACTUAL_SCALE;
	pt.y = (nodes[i]->center.y + 2*nodes[j]->center.y)/3 + 20*ACTUAL_SCALE;
	mid = new ControlPoint(pt);
	cps->push_back(mid);
    }

    cps->push_back(new ControlPoint(nodes[j]->center));
    return new VizArc(cps, this);
}

void
StructPage::OnChar( wxKeyEvent &event )
{
    if (event.m_keyCode == WXK_DELETE) {
	deleteSelectedCps();
	redraw();
	blit();
    } else {
	event.Skip();
    }
}

void
StructPage::deleteSelectedCps( void )
{
    int numNodes = nodes.size();

    for (int i = 0; i < numNodes; i++) {
	for (int j = 0; j < numNodes; j++) {
	    if (arcs[i][j]) {
		int end = arcs[i][j]->cps->size() - 1;
		for (int k = end - 1; k > 0; k--) {
		    if ((*arcs[i][j]->cps)[k]->getSelected()) {
			// delete's on it's own?
			arcs[i][j]->cps->erase(arcs[i][j]->cps->begin() + k);
			arcs[i][j]->points->DeleteNode(arcs[i][j]->points->Item(k));
		    }
		}
	    }
	}
    }
    gvpDirty = true;
    onComeForward();
}

void
StructPage::OnMouseEvent( wxMouseEvent &event )
{
    static bool dragging = false;
    static bool boxSelecting = false;
    static bool shifted = false;
    static wxPoint dragStart;
    bool screenDirty = false;
    wxPoint pt;

    event.GetPosition(&pt.x, &pt.y);
    CalcUnscrolledPosition( pt.x, pt.y, &pt.x, &pt.y );
    pt.x = (int)round(pt.x / gZoomMap[displayScale]);
    pt.y = (int)round(pt.y / gZoomMap[displayScale]);
    Selectable *pointee = itemAt(pt);

    if (event.LeftDown()) {
	shifted = event.ShiftDown(); // keep this around
	if (pointee) {
	    if (shifted) {
		pointee->toggleSelected();
		screenDirty = true;
	    }
	    else {
		dragging = true;
		dragStart.x = pt.x;
		dragStart.y = pt.y;
		if (!pointee->getSelected()) {
		    setAllSelected(false);
		    pointee->setSelected(true);
		    screenDirty = true;
		}
	    }
	} else {
	    boxSelecting = true;
	    selectBox.x = pt.x;
	    selectBox.y = pt.y;
	    selectBox.width = 0;
	    selectBox.height = 0;
	    drawSelectBox = true;
	    if (!shifted) {
		setAllSelected( false );
		screenDirty = true;
	    }
	}
    } else if (event.Dragging() && event.LeftIsDown()) {
	if (boxSelecting) {
	    selectBox.width = pt.x - selectBox.x;
	    selectBox.height = pt.y - selectBox.y;
	} else if (dragging) {
	    moveSelected( pt.x - dragStart.x, pt.y - dragStart.y );
	    dragStart.x = pt.x;
	    dragStart.y = pt.y;
	}
	screenDirty = true;
    } else if (event.LeftUp()) {
	if (boxSelecting) {
	    selectBox.width = pt.x - selectBox.x;
	    selectBox.height = pt.y - selectBox.y;
	    if (!shifted) {
		setAllSelected(false);
	    }
	    toggleSelectedInRect(selectBox);
	    screenDirty = true;
	}
	boxSelecting = shifted = dragging = drawSelectBox = false;
    } else if (event.RightDown()) {
	ControlPoint *cp = dynamic_cast<ControlPoint *>(pointee);
	if (cp) {
	    int index = -1;
	    VizArc *arc = findArcOwning(cp, index);
	    if (arc && index >= 0) {
		wxPoint where(pt.x - NEW_CP_OFFSET, pt.y - NEW_CP_OFFSET);
		arc->cps->insert( arc->cps->begin() + index,
				  new ControlPoint(where) );
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
		    wxString msg;
		    msg.sprintf("(pt.x, pt.y) = (%d, %d)\n"
				"(x0, y0)     = (%d, %d)\n"
				"(x1, y1)     = (%d, %d)\n"
				"(x, y)       = (%f, %f)\n"
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
			 fabs(y - pt.y) <= epsilon && fabs(x - pt.x) <= delta ) {
			wxPoint where(pt.x-NEW_CP_OFFSET, pt.y-NEW_CP_OFFSET);
			arc->cps->insert( arc->cps->begin() + 1,
					  new ControlPoint(where) );
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
    }

    if ( screenDirty ) {
	redraw();
	blit();
    }
    
    // Anything else to be done?
    event.Skip();
}

VizArc*
StructPage::findArcOwning( ControlPoint *cp, int& index )
{
    if (cp->arc) {
	int end = cp->arc->cps->size() - 1;
	for (int k = 1; k < end; k++) {
	    if ( (*cp->arc->cps)[k] == cp ) {
		index = k;
		return cp->arc;
	    }
	}
    }
    index = -1;
    return NULL; // nothing found
}

void
StructPage::blit( void )
{
    wxClientDC dc( this );
    PrepareDC( dc );
    blit( dc );
}

void
StructPage::blit( wxDC& dc )
{
    wxBufferedDC bdc( &dc, *content );
}

bool
StructPage::crossesNode( wxCoord x0, wxCoord x1 )
{
    int numNodes = nodes.size();

    for (int i = 0; i < numNodes; i++) {
	if ( ( ( x0 < nodes[i]->center.x-NODE_RADIUS !=
		 x1 < nodes[i]->center.x-NODE_RADIUS ) ||
	       ( x0 > nodes[i]->center.x+NODE_RADIUS !=
		 x1 > nodes[i]->center.x+NODE_RADIUS ) ) ||
	     ( ( x0 < nodes[i]->nametag.pos.x
		 + nodes[i]->nametag.size.x/2 - NODE_RADIUS !=
		 x1 < nodes[i]->nametag.pos.x - NODE_RADIUS ) ||
	       ( x0 > nodes[i]->nametag.pos.x
		 + nodes[i]->nametag.size.x/2 + NODE_RADIUS !=
		 x1 > nodes[i]->nametag.pos.x
		 + nodes[i]->nametag.size.x/2 + NODE_RADIUS ) ) )
	    return true;
    }
    return false;
}

bool
StructPage::crossesFrameEnd( wxCoord x0, wxCoord x1 )
{
    for (unsigned int i = 0; i < frameEnds.size(); i++) {
	if ((x0<frameEnds[i]->x-NODE_RADIUS!=x1<frameEnds[i]->x-NODE_RADIUS) ||
	    (x0>frameEnds[i]->x+NODE_RADIUS!=x1>frameEnds[i]->x+NODE_RADIUS))
	    return true;
    }
    return false;
}

bool
StructPage::inBounds( wxCoord x, wxCoord y )
{
    return ( x > 0 && x < getWidth() &&
	     y > 0 && y < getHeight() );
}

void
StructPage::moveSelected( int dx, int dy )
{
    int numNodes = nodes.size();

    // frame separators
    for (unsigned int i = 0; i < frameEnds.size(); i++) {
	if ( frameEnds[i]->getSelected() && inBounds(frameEnds[i]->x + dx, 1)
	     && !crossesNode(frameEnds[i]->x, frameEnds[i]->x + dx) ) {
	    frameEnds[i]->x += dx;
	}
    }
    // node name tags
    for (int i = 0; i < numNodes; i++) {
	if ( nodeNameTags[i]->getSelected()
	     && inBounds( nodeNameTags[i]->pos.x
			  + nodeNameTags[i]->size.x/2 + dx,
			  nodeNameTags[i]->pos.y
			  + nodeNameTags[i]->size.y/2 + dy ) ) {
	    if ( !crossesFrameEnd( nodeNameTags[i]->pos.x
				   + nodeNameTags[i]->size.x/2,
				   nodeNameTags[i]->pos.x + dx
				   + nodeNameTags[i]->size.x/2 ) )
		nodeNameTags[i]->pos.x += dx;
	    nodeNameTags[i]->pos.y += dy;
	}
    }
    // nodes
    for (int i = 0; i < numNodes; i++) {
	if ( nodes[i]->getSelected()
	     && inBounds(nodes[i]->center.x + dx, nodes[i]->center.y + dy) ) {
	    if (!crossesFrameEnd(nodes[i]->center.x, nodes[i]->center.x + dx))
		nodes[i]->center.x += dx;
	    nodes[i]->center.y += dy;
	    nodes[i]->tipWin->Move( nodes[i]->center.x - NODE_RADIUS,
				   nodes[i]->center.y - NODE_RADIUS );
	}
    }
    // then arcs
    for (int i = 0; i < numNodes; i++) {
	for (int j = 0; j < numNodes; j++) {
	    if (arcs[i][j]) {
		int end = arcs[i][j]->cps->size();
		assert(end > 1);
		for ( int k = 0; k < end; k++ ) {
		    if ( (*arcs[i][j]->cps)[k]->getSelected()
			 && inBounds( (*arcs[i][j]->cps)[k]->pos.x+dx,
				      (*arcs[i][j]->cps)[k]->pos.y+dy ) ) {
			if ((k && !(k == end-1)) ||
			    !crossesFrameEnd((*arcs[i][j]->cps)[k]->pos.x,
					     (*arcs[i][j]->cps)[k]->pos.x+dx))
			    (*arcs[i][j]->cps)[k]->pos.x += dx;
			(*arcs[i][j]->cps)[k]->pos.y += dy;
		    }
		}
	    }
	}
    }
    gvpDirty = true;
    onComeForward();
}

Selectable*
StructPage::itemAt( const wxPoint& pt )
{
    // XXX: This may be the best way, but is it?
    int numNodes = nodes.size();

    if (drawNodeNames) {
	// node name tags (only if they're drawn)
	for (int i = 0; i < numNodes; i++) {
	    if (nodeNameTags[i]->onMe( pt ))
		return nodeNameTags[i];
	}
    }

    // nodes
    for (int i = 0; i < numNodes; i++) {
	if (nodes[i]->onMe( pt ))
	    return nodes[i];
    }
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
    // then frame borders
    for (unsigned int i = 0; i < frameEnds.size(); i++) {
	if (frameEnds[i]->onMe( pt ))
	    return frameEnds[i];
    }
    return NULL;
}

void
StructPage::OnPaint( wxPaintEvent &WXUNUSED(event) )
{
    wxBufferedPaintDC dc( this, *content );
    //PrepareDC( dc );

    // Nothing has changed so we can skip the actual drawing and just blit 
}

void
StructPage::redraw( void )
{
    wxMemoryDC dc;
    dc.SelectObject( *content );
    dc.SetUserScale(gZoomMap[displayScale], gZoomMap[displayScale]);
    draw(dc);
}

void
StructPage::draw( wxDC& dc )
{
    dc.BeginDrawing();
    dc.SetBackground(*wxWHITE_BRUSH);
    dc.SetBrush(*wxWHITE_BRUSH);
    dc.SetFont(*wxNORMAL_FONT);
    dc.SetPen(*wxBLACK_PEN);
    dc.SetTextBackground(*wxWHITE);
    dc.SetBackgroundMode(wxTRANSPARENT);
    dc.SetTextForeground(*wxBLACK);
    dc.Clear();

    int numNodes = nodes.size();

    if (drawFrameSeps) {
	// first draw frame/chunk separators
	for (unsigned int i = 0; i < frameEnds.size(); i++) {
	    frameEnds[i]->draw(&dc);
	}
    }

    // then draw arcs
    for (int i = 0; i < numNodes; i++) {
	for (int j = 0; j < numNodes; j++) {
	    if (arcs[i][j]) arcs[i][j]->draw(&dc);
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

    if (drawNodeNames) {
	// draw node names
	for (int i = 0; i < numNodes; i++) {
	    nodeNameTags[i]->draw(&dc);
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

    dc.EndDrawing();
}

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

	    line.sprintf("strFile=%s\n", strFile.c_str());
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
		/*line.sprintf( "nodes[%d].frame=%d\n",
			      i, nodes[i]->rvId.second );
		gvp.Write(line);
		line.sprintf( "nodes[%d].nametag.name=%s\n",
			      i, nodes[i]->nametag.name.c_str() );
		gvp.Write(line);*/
		line.sprintf( "nodes[%d].nametag.pos.x=%d\n",
			      i, nodes[i]->nametag.pos.x );
		gvp.Write(line);
		line.sprintf( "nodes[%d].nametag.pos.y=%d\n",
			      i, nodes[i]->nametag.pos.y );
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
		/*line.sprintf( "frameEnds[%d].chunkBorder=%d\n",
			      i, (frameEnds[i]->chunkBorder ? 1 : 0) );
		gvp.Write(line);*/
	    }

	    gvpDirty = false;
	} else {
	    wxLogMessage(wxT("Failed to open position file for writing"));
	}
    } else {
	SaveAs();
    }
    onComeForward();
}

void
StructPage::SaveAs( void )
{
    wxString defaultGvp;
    defaultGvp.sprintf("%s.gvp", strFile.c_str());
    wxFileDialog dlg(this, "Save position file as...",
		     defaultGvp.BeforeLast('/'), defaultGvp.AfterLast('/'),
		     "All files|*"
		     "|Position Files (*.gvp)|*.gvp"
		     "|Text Files|*.txt;*.text",
		     wxSAVE | wxCHANGE_DIR | wxOVERWRITE_PROMPT,
		     wxDefaultPosition);

    dlg.SetFilterIndex(1); // show only gvps by default

    // Where do they want to save it?
    if ( dlg.ShowModal() == wxID_OK ) {
	gvpFile = dlg.GetPath();
	// save to the file
	Save();
	// put the filename wherever it needs to go
	parentFrame->SetStatusText(dlg.GetFilename(), 0);
    }
}

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
	    // try saving it
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

void
StructPage::setEndpointsSelected( bool newSelected, RVInfo::rvParent rvId )
{
    int n = nameVizNodeMap[rvId];
    int totalNodes = nodes.size();

    for ( int i = 0; i < totalNodes; i++ ) {
	if (arcs[n][i]) {
	    (*arcs[n][i]->cps)[0]->setSelected(newSelected);
	}
	if (arcs[i][n]) {
	    (*arcs[i][n]->cps)[ arcs[i][n]->cps->size()-1 ]
		->setSelected(newSelected);
	}
    }
}

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
}

void
StructPage::toggleSelectedInRect( const wxRect& rect )
{
    int numNodes = nodes.size();

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

void
StructPage::toggleViewCPs( void )
{
    drawCPs = !drawCPs;
    redraw();
    blit();
}

void
StructPage::toggleViewLines( void )
{
    drawLines = !drawLines;
    redraw();
    blit();
}

void
StructPage::toggleViewSplines( void )
{
    drawSplines = !drawSplines;
    redraw();
    blit();
}

void
StructPage::toggleViewArrowHeads( void )
{
    drawArrowHeads = !drawArrowHeads;
    redraw();
    blit();
}

void
StructPage::toggleViewNodes( void )
{
    drawNodes = !drawNodes;
    redraw();
    blit();
}

void
StructPage::toggleViewDirectLines( void )
{
    drawDirectLines = !drawDirectLines;
    redraw();
    blit();
}

void
StructPage::toggleViewFrameSeps( void )
{
    drawFrameSeps = !drawFrameSeps;
    redraw();
    blit();
}

void
StructPage::toggleViewNodeNames( void )
{
    drawNodeNames = !drawNodeNames;
    redraw();
    blit();
}

void
StructPage::toggleViewToolTips( void )
{
    drawToolTips = !drawToolTips;
    wxToolTip::Enable(drawToolTips);
}

void
StructPage::onComeForward( void )
{
    parentFrame->SetStatusText(gvpDirty?wxT("*"):wxEmptyString, 1);
    parentFrame->SetStatusText(gvpFile.length()?gvpFile:wxT("UNTITLED (and unsaved)"), 0);
}

void
StructPage::getName( wxString& name )
{
    name = ( (gvpFile != wxEmptyString) ? gvpFile : strFile );
}

void
StructPage::setScale( int newScale )
{
    displayScale = newScale;
    if (content) delete content;
    content = new wxBitmap( (int)round(canvasWidth*gZoomMap[displayScale]),
			    (int)round(canvasHeight*gZoomMap[displayScale]) );
    SetVirtualSize( (int)round(canvasWidth*gZoomMap[displayScale]),
		    (int)round(canvasHeight*gZoomMap[displayScale]) );
    do {
	wxClientDC temp(this);
	temp.SetBackground(*wxLIGHT_GREY_BRUSH);
	temp.Clear();
    } while (false);
    redraw();
    blit();
}


NameTag::NameTag( const wxPoint& newPos, const wxString& newName )
    : pos(newPos), name(newName)
{
}

void
NameTag::draw( wxDC *dc )
{
    dc->GetTextExtent(name, &size.x, &size.y);
    if (getSelected())
	dc->DrawRectangle(pos.x, pos.y, size.x, size.y);
    dc->DrawText(name, pos);
}

bool
NameTag::onMe( const wxPoint& pt )
{
    wxRect temp(pos.x, pos.y, size.x, size.y);

    return temp.Inside(pt);
}

bool
NameTag::inRect( const wxRect& rect )
{
    return rect.Inside(pos.x + size.x/2, pos.y + size.y/2);
}

VizNode::VizNode( const wxPoint& pos, RVInfo *newRvi, StructPage *newPage )
    : nametag(pos, wxT(newRvi->name.c_str()))
{
    center.x = pos.x;
    center.y = pos.y;
    rvi = newRvi;
    rvId.first = rvi->name;
    rvId.second = rvi->frame;
    page = newPage;
    tipWin = new wxWindow( newPage, -1, wxPoint(pos.x - NODE_RADIUS,
						pos.y - NODE_RADIUS),
			   wxSize(2*NODE_RADIUS, 2*NODE_RADIUS),
			   wxTRANSPARENT_WINDOW );
    tipWin->Hide();
    tipWin->SetToolTip(wxT(rvId.first.c_str()));
}

VizNode::~VizNode( void )
{
    // don't need to destroy tipWin because it's parent will
}

void
VizNode::draw( wxDC *dc )
{
    dc->DrawCircle(center, NODE_RADIUS);
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

bool
VizNode::onMe( const wxPoint& pt )
{
    int dx = (pt.x - center.x), dy = (pt.y - center.y);
    return ( hypot(dx, dy) <= NODE_RADIUS );
}

void
VizNode::setSelected( bool newSelected ) {
    selected = newSelected;
    // the selected state of the node must be the same as all the
    // control points belonging to the node
    page->setEndpointsSelected( newSelected, rvId );
    // select the node's nametag as well
    nametag.setSelected( newSelected );
}


bool
ControlPoint::onMe( const wxPoint& pt )
{
    return pos.x-1*ACTUAL_SCALE <= pt.x && pt.x <= pos.x+1*ACTUAL_SCALE
	&& pos.y-1*ACTUAL_SCALE <= pt.y && pt.y <= pos.y+1*ACTUAL_SCALE;
}

ControlPoint::ControlPoint( const wxPoint& pt )
{
    pos.x = pt.x;
    pos.y = pt.y;
    arc = NULL;
}

VizArc::VizArc( std::vector< ControlPoint* > *newCps, StructPage *newPage )
{
    cps = newCps;
    page = newPage;
    switching = false;
    conditional = false;

    for (unsigned int i = 0; i < cps->size(); i++) {
	(*cps)[i]->arc = this;
    }
    
    int numPoints = cps->size();
    points = new wxList;
    for (int i = 0; i < numPoints; i++) {
	points->Append( (wxObject*)&(*cps)[i]->pos );
    }
}

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

void
VizArc::draw( wxDC *dc )
{
    wxPen *pen = NULL;
    if (switching && conditional)
	pen = &page->bothPen;
    else if (switching)
	pen = &page->switchingPen;
    else if (conditional)
	pen = &page->conditionalPen;

    wxPen oldPen = dc->GetPen();
    dc->SetPen(*pen);

    if ( page->getViewSplines() )
	dc->DrawSpline(points);
    if ( page->getViewLines() )
	dc->DrawLines(points);
    if ( page->getViewDirectLines() )
	dc->DrawLine((*cps)[0]->pos, (*cps)[cps->size()-1]->pos);

    dc->SetPen(oldPen);
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
	dc->SetBrush(*wxBLACK_BRUSH);
	dc->DrawPolygon( 3, arrow );
	dc->SetBrush(oldBrush);
    }

    dc->SetPen(page->controlPointPen);
    int end = points->GetCount() - 1;
    assert(end > 0);
    wxNode *node = points->GetFirst();
    node = node->GetNext();
    for ( int i = 1; i < end; i++, node = node->GetNext() ) {
	assert(node != NULL);
	wxPoint *p = (wxPoint *)node->GetData();
	if ( (*cps)[i]->getSelected() )
	    dc->DrawRectangle(p->x-(3*ACTUAL_SCALE)/2, p->y-(3*ACTUAL_SCALE)/2,
			      3*ACTUAL_SCALE, 3*ACTUAL_SCALE);
	if ( page->getViewCPs() ) {
	    dc->DrawRectangle( p->x-ACTUAL_SCALE/2, p->y-ACTUAL_SCALE/2,
			       ACTUAL_SCALE, ACTUAL_SCALE );
	}
    }
    dc->SetPen(oldPen);
}


VizSep::VizSep( wxCoord newX, StructPage *newPage, bool newChunkBorder )
{
    x = newX;
    page = newPage;
    chunkBorder = newChunkBorder;
}

void
VizSep::draw( wxDC * dc )
{
    wxPen *pen = &page->frameBorderPen;
    if (chunkBorder)
	pen = &page->chunkBorderPen;

    wxPen oldPen = dc->GetPen();
    dc->SetPen(*pen);
    dc->DrawLine( x, 0, x, page->getHeight() );
    dc->SetPen(oldPen);

    if (getSelected()) {
	dc->DrawRectangle( x-1*ACTUAL_SCALE, -1*ACTUAL_SCALE,
			   3*ACTUAL_SCALE, 3*ACTUAL_SCALE );
	dc->DrawRectangle( x-1*ACTUAL_SCALE, page->getHeight()-2*ACTUAL_SCALE,
			   3*ACTUAL_SCALE, 3*ACTUAL_SCALE );
    }
}

bool
VizSep::onMe( const wxPoint& pt )
{
    return x-1*ACTUAL_SCALE <= pt.x && pt.x <= x+1*ACTUAL_SCALE
	&& 0 <= pt.y && pt.y <= page->getHeight();
}


GmtkPrintout::GmtkPrintout(StructPage * newPage, const wxChar *title)
    : wxPrintout(title)
{
    page = newPage;
}

bool
GmtkPrintout::OnPrintPage(int page)
{
    wxDC *dc = GetDC();
    if (dc)
    {
        if (page == 1)
            DrawPageOne(dc);
	else return false;
        
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

void
GmtkPrintout::GetPageInfo( int *minPage, int *maxPage,
			   int *selPageFrom, int *selPageTo )
{
    *minPage = 1;
    *maxPage = 1;
    *selPageFrom = 1;
    *selPageTo = 1;
}

bool
GmtkPrintout::HasPage(int pageNum)
{
    return (pageNum == 1);
}

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
    
    page->draw(*dc);
    //dc->SetFont(oldFont);
}

