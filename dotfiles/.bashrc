
export PYTHONUSERBASE="/data/alice/gweelden/user_python"
export PATH="$PYTHONUSERBASE/bin:$PATH"

export ALIBUILD_WORK_DIR="/data/alice/gweelden/alice/sw"
eval "`alienv shell-helper`"

export FRAGMENTATION_DIR="/data/alice/gweelden/Fragmentation"

# Set tmpdir for ccdb access inside qsub
TMPDIR=${HOME}/.globus/

source /usr/share/bash-completion/completions/git
source /usr/share/git-core/contrib/completion/git-prompt.sh

parse_git_branch() {
  git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/(\1)/'
}

environment=""
promptend="% "
if [[ x${O2_ROOT} != x ]]; then
  environment="O2"
fi
if [[ x$(echo ${LOADEDMODULES} | grep "ninja") != x ]]; then
  if [[ x${environment} != x ]]; then
    environment="${environment} "
  fi
  environment="${environment}ninja"
fi

export PS1='\A (bash) \[\e[1;35m\]\w\[\e[38;5;99m\]$(__git_ps1)\[\e[0m\] ${promptend}'
if [[ x${environment} != x ]]; then
  environment="[${environment}]"
  promptend="> "
  export PS1='\A (bash)\e[38;5;208m\]${environment} \[\e[1;35m\]\w\[\e[38;5;99m\]$(__git_ps1)\e[38;5;208m\] ${promptend}\[\e[0m\]'
fi

source ${FRAGMENTATION_DIR}/dotfiles/.aliases
source ${FRAGMENTATION_DIR}/dotfiles/.bash_aliases

