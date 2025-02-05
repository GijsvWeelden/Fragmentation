
function setuplocal {
  export FRAGMENTATION_DIR="$HOME/cernbox/Fragmentation"
  export ALIBUILD_WORK_DIR="$HOME/alice/sw"
}
function setupremote {
  export FRAGMENTATION_DIR="/data/alice/gweelden/Fragmentation"
  export ALIBUILD_WORK_DIR="/data/alice/gweelden/alice/sw"
  export PYTHONUSERBASE="/data/alice/gweelden/user_python"
  export PATH="$PYTHONUSERBASE/bin:$PATH"
  eval "`alienv shell-helper`"
}

source ${FRAGMENTATION_DIR}/dotfiles/.promptsetup
source ${FRAGMENTATION_DIR}/dotfiles/.aliases
