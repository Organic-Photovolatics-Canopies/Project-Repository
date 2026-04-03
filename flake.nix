{
  description = "An FHS shell with conda and cuda.";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = {
    self,
    nixpkgs,
    home-manager,
  }: let
    pkgs = import nixpkgs {
      system = "x86_64-linux";
      config.allowUnfree = true;
    };

    # Conda installs it's packages and environments under this directory
    installationPath = "~/.conda";

    # Downloaded Miniconda installer
    minicondaScript = pkgs.stdenv.mkDerivation rec {
      name = "miniconda-${version}";
      version = "4.3.11";
      src = pkgs.fetchurl {
        url = "https://repo.continuum.io/miniconda/Miniconda3-${version}-Linux-x86_64.sh";
        sha256 = "1f2g8x1nh8xwcdh09xcra8vd15xqb5crjmpvmc2xza3ggg771zmr";
      };
      # Nothing to unpack.
      unpackPhase = "true";
      # Rename the file so it's easier to use. The file needs to have .sh ending
      # because the installation script does some checks based on that assumption.
      # However, don't add it under $out/bin/ becase we don't really want to use
      # it within our environment. It is called by "conda-install" defined below.
      installPhase = ''
        mkdir -p $out
        cp $src $out/miniconda.sh
      '';
      # Add executable mode here after the fixup phase so that no patching will be
      # done by nix because we want to use this miniconda installer in the FHS
      # user env.
      fixupPhase = ''
        chmod +x $out/miniconda.sh
      '';
    };

    # Wrap miniconda installer so that it is non-interactive and installs into the
    # path specified by installationPath
    conda =
      pkgs.runCommand "conda-install"
      {buildInputs = [pkgs.makeWrapper minicondaScript];}
      ''
        mkdir -p $out/bin
        makeWrapper                            \
          ${minicondaScript}/miniconda.sh      \
          $out/bin/conda-install               \
          --add-flags "-p ${installationPath}" \
          --add-flags "-b"
      '';

    vscode = pkgs.vscode-with-extensions.override {
      vscode = pkgs.vscodium.fhs;
      vscodeExtensions = with pkgs.vscode-extensions; [
        vscodevim.vim
        streetsidesoftware.code-spell-checker
        ms-toolsai.jupyter
        ms-python.python
        ms-python.debugpy
        ms-python.pylint
        detachhead.basedpyright
      ];
    };
  in {
    devShell.x86_64-linux =
      (pkgs.buildFHSEnv {
        name = "conda";
        targetPkgs = pkgs: (
          with pkgs; [
            vscode
            autoconf
            binutils
            conda
            cudatoolkit
            curl
            freeglut
            gcc13
            git
            gitRepo
            gnumake
            gnupg
            gperf
            libGLU
            libGL
            libselinux
            linuxPackages.nvidia_x11
            m4
            ncurses5
            procps
            stdenv.cc
            unzip
            util-linux
            wget
            libice
            libsm
            libx11
            libxext
            libxi
            libxmu
            libxrandr
            libxrender
            libxv
            zlib
          ]
        );
        profile = ''
          # cuda
          export CUDA_PATH=${pkgs.cudatoolkit}
          # export LD_LIBRARY_PATH=${pkgs.linuxPackages.nvidia_x11}/lib
          export EXTRA_LDFLAGS="-L/lib -L${pkgs.linuxPackages.nvidia_x11}/lib"
          export EXTRA_CCFLAGS="-I/usr/include"

          # conda
          export PATH=${installationPath}/bin:$PATH
          # Paths for gcc if compiling some C sources with pip
          export NIX_CFLAGS_COMPILE="-I${installationPath}/include"
          export NIX_CFLAGS_LINK="-L${installationPath}lib"
          # Some other required environment variables
          export FONTCONFIG_FILE=/etc/fonts/fonts.conf
          export QTCOMPOSE=${pkgs.libx11}/share/X11/locale
        '';
      }).env;
  };
}
