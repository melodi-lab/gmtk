# Copyright 1999-2014 Gentoo Foundation
# Distributed under the terms of the GNU General Public License v2
# $Header: $

EAPI=5

DESCRIPTION="The Graphical Models Toolkit"
HOMEPAGE="https://j.ee.washington.edu/trac/gmtk"
SRC_URI="http://melodi.ee.washington.edu/downloads/gmtk/${P}.tar.gz"

LICENSE="OSL-3"
SLOT="0"
KEYWORDS="~amd64 ~x86"
IUSE="hdf5 wxwidgets"

RDEPEND="
	hdf5? ( sci-libs/hdf5 )
	wxwidgets? ( >=x11-libs/wxGTK-3.0.0 )"
DEPEND="${RDEPEND}"
