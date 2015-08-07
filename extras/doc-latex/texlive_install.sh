#!/bin/bash
# By default, the main CTAN mirror is used to fetch packages from. If this
# is distant from you, or you have a faster local mirror simply override
# this by setting the CTAN_MIRROR_URL environmental variable outside this
# script.
# A list of mirrors can be found at http://ctan.org/mirrors

# e.g. before running this script do:
#   export CTAN_MIRROR_URL='http://mirror.aarnet.edu.au/pub/CTAN'

wget \
  --continue \
  --directory-prefix /tmp \
  ${CTAN_MIRROR_URL:-'http://mirror.ctan.org'}/systems/texlive/tlnet/install-tl-unx.tar.gz
tar \
  --extract \
  --gunzip \
  --directory /tmp \
  --file /tmp/install-tl-unx.tar.gz

# Install texlive using the supplied texlive.profile (this just installs a
# basic LaTeX environment
/tmp/install-tl-*/install-tl \
  -repository ${CTAN_MIRROR_URL:-'http://mirror.ctan.org'}/systems/texlive/tlnet \
  -no-gui \
  -profile texlive.profile

# Install packages required by the project
packages=(
  latex-bin
  oberdiek
  url
  graphics
  pdftex-def
  colortbl
  hyperref
  xcolor
  tools
  hanging
  changepage
  geometry
  enumitem
  latexmk
)
tlmgr \
  -repository ${CTAN_MIRROR_URL:-'http://mirror.ctan.org'}/systems/texlive/tlnet \
  install \
    ${packages[@]}

