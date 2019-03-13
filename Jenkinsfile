pipeline {

  options {
    buildDiscarder(logRotator(numToKeepStr: '5'))
  }

  environment {
        PATH="/opt/anaconda/bin:$PATH"
  }

  agent {
    node {
      label 'ci-community-docker'
    }
  }

  stages {

    stage('Package & Dockerize') {
      steps {
        withMaven( maven: 'apache-maven-3.0.5' ) {
            sh 'mvn -B deploy'
        }
      }
    }
  }
}
