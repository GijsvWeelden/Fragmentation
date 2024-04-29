export ALIBUILD_WORK_DIR="$HOME/alice/sw"
eval "`alienv shell-helper`"
# Change prompt to show environment name
# Can be made fancier by checking $LOADEDMODULES for O2/ninja or other libraries. Could also show git branch, etc.
if [ x$O2_ROOT != x ]
then
export PS1='%T %F{208}[O2 ninja] %F{magenta}%~%F{105}${vcs_info_msg_0_}%F{208} %#> %f'
else
export PS1='%T %F{magenta}%~%F{105}${vcs_info_msg_0_}%f %# '
fi
export CLICOLOR=1
export LSCOLORS=gxgxBxDxCxEgEdxbxgxcxd
# Version control settings
autoload -Uz compinit && compinit
autoload -Uz vcs_info
precmd() { vcs_info }
zstyle ':vcs_info:git*' formats ' (%b)'
setopt PROMPT_SUBST
# Aliases
alias o2="alienv enter --shellrc O2Physics/latest ninja/latest"
alias unfold="alienv enter --shellrc O2Physics/latest ninja/latest RooUnfold/latest"
alias ls="ls -GF"
alias ll="ls -GF -lht"
alias lg=LG
alias lgalt=LGalt
alias nikhef="ssh -A gweelden@login.nikhef.nl"
alias stoom="ssh -A stbc-i1"
alias stoomboot="ssh -A stbc-i2"
alias stoomboter="ssh -A stbc-i3"
alias dep="$O2PHYSICS_ROOT/share/scripts/find_dependencies.py -t" # Finds table producer
alias applyformatting="clang-format -style=file -i" # Usage: format <file> # Applies O2Physics formatting
alias copytrainresults="~/Documents/Fragmentation/scripts/copyTrainResults" # Requires to be in O2 environment

# ls show if dir is git repo
LGalt() {
  paste <(CLICOLOR_FORCE=true ls -ld *) <(for i in *; do if [ -d "$i"/.git ] ; then echo "($(git --git-dir="$i"/.git symbolic-ref --short HEAD))"; else echo; fi;
  done)
}
LG() {
  if [ x$1 != x ]; then
  paste <(CLICOLOR_FORCE=true ls -ld $1) <(for i in $1; do if [ -d "$i"/.git ] ; then echo "($(git --git-dir="$i"/.git symbolic-ref --short HEAD))"; else echo; fi; done)
  else
  LGalt
  fi
}

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

