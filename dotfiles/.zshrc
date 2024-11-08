
# export ALIBUILD_WORK_DIR, FRAGMENTATION_DIR wherever you call this from
# eval "`alienv shell-helper`" wherever you call this from
# equip with conda loader if needed

# Change prompt to show environment name etc.
hostname=${HOST:0:7}
hostname="${hostname} "
if [[ $SSH_CONNECTION == "" ]]; then
  hostname=""
fi

environment=""
promptend="%f%% "
if [[ x${O2PHYSICS_ROOT} != x ]]; then
  environment=${${O2PHYSICS_ROOT##*/}%-*}
  environment="%F{208}[${environment}]"

  promptend="%F{208}%#> %f"
fi

gitstring="%F{105}${vcs_info_msg_0_}"
workdir="%F{magenta}%~"

export PS1='%T ${hostname}${environment}${workdir} ${gitstring} ${promptend}'
# export PS1='%T ${hostname}%F{208}${environment}%F{magenta}%~%F{105}${vcs_info_msg_0_} ${promptend}'

export CLICOLOR=1
export LSCOLORS=gxgxBxDxCxEgEdxbxgxcxd

# Version control settings
autoload -Uz compinit && compinit
autoload -Uz vcs_info
precmd() { vcs_info }
zstyle ':vcs_info:git*' formats ' (%b)'
setopt PROMPT_SUBST

source ${FRAGMENTATION_DIR}/dotfiles/.aliases
source ${FRAGMENTATION_DIR}/dotfiles/.zsh_aliases
