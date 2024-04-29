export ALIBUILD_WORK_DIR="$HOME/alice/sw"
export FRAGMENTATION_DIR="$HOME/cernbox/Fragmentation"
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

source .aliases
source .zsh_aliases

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

