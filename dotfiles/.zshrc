
function setuplocal {
  export FRAGMENTATION_DIR="$HOME/cernbox/Fragmentation"
  export ALIBUILD_WORK_DIR="$HOME/alice/sw"
  eval "`alienv shell-helper`"
  # >>> conda initialize >>>
  # !! Contents within this block are managed by 'conda init' !!
  __conda_setup="$('/Users/gijsvanweelden/anaconda3/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
  if [ $? -eq 0 ]; then
      eval "$__conda_setup"
  else
      if [ -f "/Users/gijsvanweelden/anaconda3/etc/profile.d/conda.sh" ]; then
          . "/Users/gijsvanweelden/anaconda3/etc/profile.d/conda.sh"
      else
          export PATH="/Users/gijsvanweelden/anaconda3/bin:$PATH"
      fi
  fi
  unset __conda_setup
  # <<< conda initialize <<<
}

function setupremote {
  export FRAGMENTATION_DIR="/data/alice/gweelden/Fragmentation"
  export ALIBUILD_WORK_DIR="/data/alice/gweelden/alice/sw"
  export LS_COLORS='di=36:ln=1;34:so=1;31:pi=1;33:ex=1;32:bd=1;34;46:cd=1;34;43:su=0;41:sg=0;46:tw=0;42:ow=0;43'
  source /cvmfs/alice.cern.ch/etc/login.sh
  # Load Git completion
  zstyle ':completion:*:*:git:*' script ~/.zsh/git-completion.bash
  fpath=(~/.zsh $fpath)
  autoload -Uz compinit && compinit

}

if [[ $SSH_CONNECTION == "" ]]
then setuplocal
else setupremote
fi

source $FRAGMENTATION_DIR/dotfiles/.setup
source $FRAGMENTATION_DIR/dotfiles/.aliases
